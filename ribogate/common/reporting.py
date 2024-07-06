import numpy as np
import math
import copy
import xlsxwriter
from io import BytesIO

from .utility import arg_sort
import config

def sort_individuals_by_fitness_name(individuals, fitness_name, is_order_ascending):
    vals = [ind.potential_fitnesses[fitness_name] for ind in individuals]
    sorted_idxs = arg_sort(vals)
    if not is_order_ascending:
        sorted_idxs = sorted_idxs[::-1]
    sorted_inds = [individuals[idx] for idx in sorted_idxs]
    
    return sorted_inds

def rank_by_marginal_diversity(distance_matrix, num_inds):
    '''
    Ranks the individuals in such a way that the Nth individual is the one that is the furthest away (in phenotype space) from the previous N-1 individuals. Uses a greedy algorithm. 
    
    Parameters:
        num_inds (int): The number of individuals being ranked
        distance_matrix (numpy float array of size num_inds X num_inds): Entry (i,j) the phenotypic distance between individuals i and j
    
    Returns:
        ranks (list of ints): The ith entry in this list stores the rank of the ith individual in the population on which the distance_matrix was calculated    
    '''

    if num_inds == 0: # no individuals
        return None 
    elif num_inds == 1:
        ranks = [0] # If there is only one individual, theh ranking it is trivial
    else:        
        sorted_idxs = list() # The ith entry of this list stores the index of the individual with rank i   
        
        # We start by finding the two furthest individuals in the distance matrix and placing them in the sorted_idxs list
        furthest_pair_idx = np.argmax(distance_matrix) # This is the scalar index of the pair in a 2d matrix
        sorted_idxs += [furthest_pair_idx % num_inds, furthest_pair_idx // num_inds] # From this scalar pair index, we get the index of each element of the pair
        
        # We then loop over all unsorted individuals and place the one that is furthest to all sorted individuals
        while len(sorted_idxs) < num_inds:
            candidate_distances = list()
            unsorted_idxs = [idx for idx in range(num_inds) if idx not in sorted_idxs]
            
            # For each unsorted individual, retrieve the distance to the closest sorted individual. This is the 'candidate_distance'
            for unsorted_idx in unsorted_idxs:
                marginal_distances = list()
                for sorted_idx in sorted_idxs:
                    marginal_distances.append(distance_matrix[unsorted_idx, sorted_idx])
                candidate_distances.append( min(marginal_distances)  )
                
            # Then we find the unsorted individual point that is furthest to all sorted individuals (i.e. the one with the highest candidate_distance value)                
            furthest_candidate_idx = arg_sort(candidate_distances)[-1]
            sorted_idxs.append( unsorted_idxs[furthest_candidate_idx] )

        ranks = [None]*num_inds
        for rank, idx in enumerate(sorted_idxs):
            ranks[idx] = rank
    return ranks


def determine_reported_individuals_and_marginal_diversity_rank(individuals, criteria, distance_matrix_func, distance_matrix_func_params):
            
    # Determine which individuals are acceptable to report to the user
    reported_individuals = list()
    for ind in individuals:
        is_reported = True
        for criterion in criteria:
            if ind.potential_fitnesses[criterion['name']] < criterion['floor'] or ind.potential_fitnesses[criterion['name']] > criterion['ceil']:
                is_reported = False
                break
        if is_reported:
            reported_individuals.append(ind)      
        ind.potential_fitnesses['is_reported'] = is_reported

    # Rank the reported individuals by marginal diversity
    distance_matrix = distance_matrix_func(reported_individuals, distance_matrix_func_params) 
    ranks = rank_by_marginal_diversity(distance_matrix, len(reported_individuals))
    for i in range(len(reported_individuals)):
        reported_individuals[i].potential_fitnesses['marginal_diversity_rank'] = -ranks[i] # Negative value so that worse ranking individuals have lower value
            
    # Individuals that were not deemed acceptable are not actively ranked and are instead given the worst possible rank        
    for ind in individuals:
        if not ind.potential_fitnesses['is_reported']:               
            ind.potential_fitnesses['marginal_diversity_rank'] = -len(reported_individuals) # This is the worst possible md rank
                    
    return individuals


def generate_srz_dna_construct(ribozyme_sequence, construct_segment_names, construct_segments, environment):
    '''
    Packages the ribozyme sequence inside of a construct containing other possible segments. 
    It is important to note that the construct segments are distinct from the ribozyme segments such as S1A, C1, etc.
    In the construct, the ribozyme is a single segment. There may be other segments such as CTE and tRNA that help the ribozyme function in a biological environment
    Note however, that these additional segements play no role in the EA. They are just added after the fact.

    Parameters:
        ribozyme_sequence (string)
        construct_segments_names (dict of list of strings): For each environment type, specifies the names of segments that appear in the construct (ordered from 5' to 3')
        construct_segments (dict of dicts). Specifies the name, sequence, and description for each possible construct segment
        environment (string): The environment that the user selected in the request page. 

    Returns:
        ribozyme_construct (list of dicts). Each dict contains the sequence and color of a construct segment

    '''

    '''
    Most of the construct segment sequences are known beforehand and constant for all individuals
    However, the ribozyme sequence depends on the individual
    Also, if the construct contains a T7 promoter, we need to ensure that the transcription start site should have at least 2 Gs
    '''
    construct_segments_complete = copy.deepcopy(construct_segments)
    construct_segments_complete['selective ribozyme']['sequence'] = str(ribozyme_sequence)

    promoter_sequence = ''
    if environment == 'in_vitro_default_T7':
        promoter_sequence = construct_segments_complete['promoter']
        # Transcription start site should have at least 2 Gs
        if ribozyme_sequence[0] != 'G': 
            promoter_sequence += 'GG'
        elif ribozyme_sequence[1] != 'G':
            promoter_sequence += 'G'
    elif 'in_vitro_custom_promoter' in environment:
        promoter_sequence = environment.split('_')[-1]
    construct_segments_complete['promoter']['sequence'] = promoter_sequence

    ribozyme_construct = []
    for name in construct_segment_names[environment]:
        # The 'text' and 'color' keys are expected by the job_details.html tmeplate
        ribozyme_construct.append({'text': construct_segments_complete[name]['sequence'], 'color':  construct_segments_complete[name]['color']}) 

    return ribozyme_construct

def generate_srz_dna_construct_header(construct_segment_names, construct_segments, environment):
    '''
    Packages the construct names and colors in order to form a color legend that will be eventually displayed on the client's browser and downloaded spreadsheet
    '''

    if len(construct_segment_names[environment]) == 1:
        header_start_text = "Selective ribozyme sequence "
    else:
        header_start_text  = "Selective ribozyme construct sequence ( "
    
    header = [{'text': header_start_text , 'color': "black"}]
    if len(construct_segment_names[environment]) > 1:
        segment_seperation_char = '-' 
        for i, name in enumerate(construct_segment_names[environment]):
            if i == len(construct_segment_names[environment]) - 1:
                # We don't put a separation character after the last segment
                segment_seperation_char = '' 
                
            # The 'text', 'color' and 'description' keys are expected by the job_details.html tmeplate
            header.append({'text': name + segment_seperation_char, 'color': construct_segments[name]['color'], 'description': construct_segments[name]['description']}) 
        header_end_text = ")"
        header.append( {'text': header_end_text , 'color': "black"} )
    return header

# I don't like putting this function here, but I prefer that to importing it from TS.folding_preparation
def get_input_states(num_inputs):
    '''
    Returns the input columns of the truth table for a given number of inputs
    Currently harcoded for 1-3 inputs.
    If we want to handle an arbitrary number of inputs we will have to generate the input_states dynamically
    '''

    if num_inputs == 1:
        input_states = [(0, ), (1, )]
    elif num_inputs == 2:
        input_states = [(0, 0), (0, 1), (1,0), (1,1)]
    elif num_inputs == 3:
        input_states = [(0, 0, 0), (0, 0, 1), (0, 1,0), (0, 1,1), (1, 0, 0), (1, 0, 1), (1, 1,0), (1, 1,1)]    

    return input_states

def get_state_names(sub_app_name, num_states):
    # Prepare the state names that are shown in the truth table    
    if sub_app_name == 'ts':        
        num_inputs = int(math.log2(num_states))
        input_states = get_input_states(num_inputs)
        state_names = [[str(bit) for bit in input_state] for input_state in input_states]
    elif sub_app_name == 'tc':
        if num_states <= 2:
            raise ValueError("In TriCleaver, there must be at least 3 states. Something went wrong upstream.")
        state_names = ["SRz "]
        for i in range(num_states-2):
            state_names.append("SRz + WT configuration " + str(i) + " ")
        state_names.append("SRz + mutant ")  
        
    return state_names

def preformat_segments(segments, attribute_to_display, sub_app_name, should_add_ampersand=False, should_add_dash=False):
    '''
    Wrap relevant segment properties in a list of dictionaries so that they may be JSON serialized
    This list of dictionaries is also used to write rich spreadsheet text
    '''
    if sub_app_name == 'tc':
        num_folded_strands = 2
    elif sub_app_name == 'ts':
        num_folded_strands = 1
        
    preformatted_segments = list()
    for i in range(num_folded_strands): # Don't show the non-folded strands (i.e inputs for TS and repeats for TC)
        strand_segments = [seg for seg in segments if seg.strand == i]
        for j, seg in enumerate(strand_segments):
            if attribute_to_display == 'name':
                text_attribute = seg.name
            elif attribute_to_display == 'verbose_name':
                text_attribute = config.verbose_seg_names[seg.name]
            elif attribute_to_display == 'sequence':
                text_attribute = seg.seq

            preformatted_seg = {'text': text_attribute, 'color': config.seg_colors[seg.name]}
            preformatted_segments.append(preformatted_seg)
            
            # Currently, used only for TC and TS spreadsheet
            if should_add_dash and (i != num_folded_strands-1 or j != len(strand_segments)-1):
                preformatted_dash = {'text': '-', 'color': 'black'}
                preformatted_segments.append(preformatted_dash)
                
        # Currently, used only for TC spreadsheet
        if should_add_ampersand and i != num_folded_strands-1:
            preformatted_ampersand = {'text': '&', 'color': 'black'}
            preformatted_segments.append(preformatted_ampersand)
            
    return preformatted_segments

def get_marked_cut_region(extended_rzbs, segments):
    '''
    inputs:
        extended_rzbs: The RzBS sequence, along with the upstream_rzbs and the dowstream_rzbs regions.
    '''
   
    # The len upstream_rzbs_seq may vary between individuals with rzbs_position near 5' end of downstream if we used the POST rzbs_extension_len 
    upstream_rzbs_seq = [seg for seg in segments if seg.name == 'UPSTREAM_RZBS'][0].seq
    s3b_seq = [seg for seg in segments if seg.name == 'S3B'][0].seq
    cut_position = len(upstream_rzbs_seq) + len(s3b_seq) + 1

    marked_cut_region = extended_rzbs[0:cut_position] + "|" + extended_rzbs[cut_position:]

    return marked_cut_region


def write_rich_text_spreadsheet(workbook, worksheet, cell, formatted_segments):
    
    formatted_string = list()
    for segment in formatted_segments:            
        segment_format = workbook.add_format({'font_color': segment['color'], 'font_name':config.spreadsheet_font_name})
        segment_sequence = segment['text']
        formatted_string = formatted_string + [segment_format, segment_sequence]   
    if len(formatted_segments)> 1:
        worksheet.write_rich_string(cell, *formatted_string)    
    else:
        worksheet.write_string(cell, segment_sequence, segment_format) # write_rich_string fails if there is only one formatted segment   

def export_spreadsheet_bytes(ranked_inds, job_info, sub_app_name):
    output = BytesIO()
    workbook = xlsxwriter.Workbook(output)
    
    # Sheet 1: Sequences
    worksheet = workbook.add_worksheet("Sequences")
    worksheet.write_string('A1', "Rank") # Maybe: Change 'Rank' to 'Recommended testing order'

    if sub_app_name == 'tc':
        # Write header row
        construct_header = generate_srz_dna_construct_header(config.construct_segment_names, config.construct_segments, job_info['ENVIRONMENT'])    
        write_rich_text_spreadsheet(workbook, worksheet, 'B1', construct_header)
        worksheet.write_string('C1', "Cut region")
        worksheet.write_string('D1', "Activity score")
        worksheet.write_string('E1', "Inativity score")
        worksheet.write_string('F1', "Mismatch score")
 
        # Write and color the sequence of the entire of packaged ribozyme
        for rank, ind in enumerate(ranked_inds):             
            ribozyme_sequence = str(ind.strands[0].replace('U','T'))
            ribozyme_construct = generate_srz_dna_construct(ribozyme_sequence, config.construct_segment_names, config.construct_segments, job_info['ENVIRONMENT'])   
    
            worksheet.write_string('A' + str(rank+2), str(rank+1) )
            write_rich_text_spreadsheet(workbook, worksheet, 'B' + str(rank+2), ribozyme_construct)
            
            marked_cut_region = get_marked_cut_region(ind.strands[1], ind.segments)
            marked_cut_region = marked_cut_region.replace('U','T')
            worksheet.write_string('C' + str(rank+2), marked_cut_region)
            
            on_score, off_score, mismatch_score = [format(ind.potential_fitnesses[name], f'.{3}g') for name in ['ON', 'OFF', 'mismatch']] # 3 sig figs    
            worksheet.write_string('D' + str(rank+2), on_score)
            worksheet.write_string('E' + str(rank+2), off_score)
            worksheet.write_string('F' + str(rank+2), mismatch_score)

    elif sub_app_name == 'ts':
        # Write header rpw
        worksheet.write_string('B1', "Ribogate sequence" )
        # The number of inputs can very between 1 and 3 depending on the job, so we can't hardcode column letters
        column_letter_code = ord('C')
        for i in range(job_info['NUM_INPUTS']):
            worksheet.write_string(chr(column_letter_code + i) + '1', "Input " + str(i+1) + " sequence" )
        
        worksheet.write_string(chr(column_letter_code + job_info['NUM_INPUTS'])  + '1', "Activity score")
        worksheet.write_string(chr(column_letter_code + job_info['NUM_INPUTS'] + 1) + '1', "Inativity score")

        for rank, ind in enumerate(ranked_inds):     
            worksheet.write_string('A' + str(rank+2), str(rank+1) )

            column_letter_code = ord('B')
            for i in range(job_info['NUM_INPUTS']+1):
                worksheet.write_string(chr(column_letter_code + i) + str(rank+2), ind.strands[i])
                
            on_score, off_score = [format(ind.potential_fitnesses[name], f'.{3}g') for name in ['ON', 'OFF']]
            worksheet.write_string(chr(column_letter_code + job_info['NUM_INPUTS'] + 1) + str(rank+2) , on_score)
            worksheet.write_string(chr(column_letter_code + job_info['NUM_INPUTS'] + 2) +  str(rank+2), off_score)
                


    # Sheet 2: Predicted structures    
    worksheet2 = workbook.add_worksheet("Predicted structures") 
    font_format = workbook.add_format({'font_name':config.spreadsheet_font_name})
    
    worksheet2.write_string('A1', "Segment color legend")
    preformatted_segments = preformat_segments(ind.segments, 'verbose_name', sub_app_name, False, True)
    write_rich_text_spreadsheet(workbook, worksheet2, 'B1', preformatted_segments)

    delimeter = "###"
    row_idx = 3
    for rank, ind in enumerate(ranked_inds):

        if sub_app_name == 'ts':
            sequence = ind.strands[0]
            constraints = [task[1] for task in ind.tasks['logic']]
            preformatted_segments = segments_js = preformat_segments(ind.segments, 'sequence', 'ts')
            structures = [structure for structure in ind.mfe_dot_bracket_structures]
        elif sub_app_name == 'tc':
            sequence = ind.strands[0] + "&" + ind.strands[1]
            constraints = [task[1] for task in ind.tasks]
            preformatted_segments = preformat_segments(ind.segments, 'sequence', 'tc', True)
            
            structures = [structure for structure in ind.mfe_dot_bracket_structures]
            if config.should_report_unhybridized_srz:   
                structures[0] = ind.mfe_db_unhybridized_srz 
            
        worksheet2.write_string('A' + str(row_idx), delimeter + " BEGIN INDIVIDUAL WITH RANK " + str(rank) + ' ' + delimeter )
        worksheet2.write_string('B' + str(row_idx), '#'*len(sequence), font_format )
        row_idx += 1            

        num_states = len(ind.mfe_dot_bracket_structures)
        state_names = get_state_names(sub_app_name, num_states)       
        for i in range(num_states):
            state_name = ''.join(state_names[i]) # The join method handles the list of bits of TS and doesn't effect the existing string for TC
            worksheet2.write_string('A' + str(row_idx), state_name, font_format)
            row_idx += 1
            
            worksheet2.write_string('A' + str(row_idx), "Secondary structure", font_format)
            worksheet2.write_string('B' + str(row_idx), ind.mfe_dot_bracket_structures[i], font_format)
            row_idx += 1
            
            worksheet2.write_string('A' + str(row_idx), "Constraint", font_format)
            worksheet2.write_string('B' + str(row_idx), constraints[i], font_format)
            row_idx += 1            

            worksheet2.write_string('A' + str(row_idx), "Sequence", font_format)
            #worksheet2.write_string('B' + str(row_idx), sequence)
            write_rich_text_spreadsheet(workbook, worksheet2, 'B' + str(row_idx), preformatted_segments)
            row_idx += 2
            
        worksheet2.write_string('A' + str(row_idx), delimeter + " END INDIVIDUAL WITH RANK " + str(rank) + ' ' + delimeter , font_format)
        worksheet2.write_string('B' + str(row_idx), '#'*len(sequence) , font_format)
        row_idx += 2
        
    workbook.close()
    output.seek(0)
    
    return output 

def export_spreadsheet(individuals, job_info, sub_app_name, filename):
    
    reportable_inds = [ind for ind in individuals if ind.potential_fitnesses['is_reported']]
    ranked_inds = sort_individuals_by_fitness_name(reportable_inds, 'marginal_diversity_rank', False)

    output = export_spreadsheet_bytes(ranked_inds, job_info, sub_app_name)
    with open(filename, 'wb') as f:
        f.write(output.read())
        
