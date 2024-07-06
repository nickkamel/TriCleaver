from .tc.substrate_processing import extract_every_rzbs, count_tandem_repeats
import config
from .common.utility import convert_to_float
from .tc.ribozyme_templates import RibozymeTemplate
   
def validate_wt_mutant_seqs(wt_seq, mutant_seq):  
    '''
    Validate the wild-type and mutant sequences
    '''
        
    errors = {'wt':[], 'mutant':[] }

    if not wt_seq:
        errors['wt'].append('Please specify the sequence of the wild-type strand')
        return None, errors
    
    wt_seq = wt_seq.replace(' ','').replace('\n','').replace('a','A').replace('g','G').replace('c','C').replace('t','U').replace('T','U')

    if not mutant_seq:
        errors['mutant'].append('Please specify the sequence of the mutant strand')
        return None, errors
    
    mutant_seq = mutant_seq.replace(' ','').replace('\n','').replace('a','A').replace('g','G').replace('c','C').replace('t','U').replace('T','U')

    if wt_seq.count('-') != 2:
        errors['wt'].append("The repeats region must be delineated by one dash on each side. e.g. AAAA-CAGCAGCAG-AAAA")
        return None, errors 

    upstream_wt, repeats_wt, downstream_wt = wt_seq.split('-')

    if mutant_seq.count('-') != 2:
        errors['mutant'].append("The repeats region must be delineated by one dash on each side. e.g. AAAA-CAGCAGCAG-AAAA")
        return None, errors
    
    upstream_mutant, repeats_mutant, downstream_mutant = mutant_seq.split('-')

    valid_chars = set("AUGTCaugtc-")
    if not set(wt_seq).issubset(valid_chars):
        errors['wt'].append("Wild-type sequence contains invalid characters. Allowed characters are A, U, G, T, C, as well a dash (-) on each side of the repeats segment")
        return None, errors
    if not set(mutant_seq).issubset(valid_chars):
        errors['mutant'].append("Mutant sequence contains invalid characters. Allowed characters are A, U, G, T, C, as well a dash (-) on each side of the repeats segment")
        return None, errors

    if upstream_wt != upstream_mutant:
        error_message = "The sequence of the region upstream of the repeats must be identical between the wild-type and mutant."
        errors['mutant'].append(error_message)
        errors['wt'].append(error_message)
        return None, errors
    
    upstream_seq = upstream_wt

    if downstream_wt != downstream_mutant:
        error_message = "The sequence of the region downstream of the repeats must be identical between the wild-type and mutant."
        errors['mutant'].append(error_message)
        errors['wt'].append(error_message)
        return None, errors
    
    downstream_seq = downstream_wt

    if (len(mutant_seq) - 3) < len(wt_seq): # This '3' will have to be replaced if we want to deal with repeats of size 4,5,6
        errors['mutant'].append("The mutant must contain at least one more repeat than the wild-type.")
        errors['wt'].append("The wild-type must contain at least one less repeat than the mutant.")
        return None, errors

    data = {'wt_seq':wt_seq, 'mutant_seq':mutant_seq, 'downstream_seq':downstream_seq, 'upstream_seq':upstream_seq, 'repeats_wt':repeats_wt, 'repeats_mutant':repeats_mutant}
    return data, None

def validate_repeats(repeats_wt, repeats_mutant, repeat_unit_len):
    '''
    Determine the # of wt and mutant repeats
    '''    
    errors = {'wt':[], 'mutant':[] }

    repeat_unit_wt, num_repeats_wt = count_tandem_repeats(repeats_wt, repeat_unit_len)
    repeat_unit_mutant, num_repeats_mutant = count_tandem_repeats(repeats_mutant, repeat_unit_len)
    if repeat_unit_wt != repeat_unit_mutant:
        error_message = "The type of repeat (e.g. GCG, CAG) must be the same between the wild-type and mutant."
        errors['mutant'].append(error_message)
        errors['wt'].append(error_message)
        return None, errors
    if num_repeats_wt >= num_repeats_mutant:
        errors['mutant'].append("The number of mutant repeats must exceed the number of wild-type repeats.")
        return None, errors
    if num_repeats_wt < 5:
        errors['wt'].append("The wild type must have at least 5 repeats.")
        return None, errors
    repeat_unit = repeat_unit_wt

    '''
    In some cases (e.g. HTT, MJD) there are a lot of mutant repeats relative to the wild-type. 
    In one sense this is good: the higher the ratio between them, the easier it is for the SENSOR to differentiate between them. 
    However, at some point there is overkill. The SENSOR will be very long which will grind folding to a halt without further increasing selectivitly. 
    Therefore, we cap the SENSOR length.
    '''
    if num_repeats_mutant - num_repeats_wt > config.template_params_tc['MAX_MUTANT_WT_DIFFERENCE']:
        num_repeats_mutant = num_repeats_wt + config.template_params_tc['MAX_MUTANT_WT_DIFFERENCE']
        
    data = {'num_repeats_wt':num_repeats_wt, 'num_repeats_mutant':num_repeats_mutant, 'repeat_unit':repeat_unit}
    return data, None

def validate_potential_cut_sites(ribozyme_type, downstream_seq):
    '''
    Extract potential ribozyme binding sites (rzbs) and their associated cut site motifs (NUH)
    '''

    errors = {'cut_site_motifs':[]}

    '''
    The potential ribozyme binding sites can vary depending on which template is used
    Currently, this only occurs because the extended hammerhead is a bit longer so RzBSs near the end of the downstream seq might not fit in
    We therefore instantiate the selected ribozyme type's template in order to access its stem 3 length and its RzBS length
    This template instance is not used anywhere else in the code (we instantiate it again later in the TriCleaver EA)
    Note we just specify a blank repeats sequence since that doesn't effect any of the parameters we need to query here
    '''
    template = RibozymeTemplate(ribozyme_type, config.template_params_tc, downstream_seq, '')    

    rzbs_positions, rzbs_cut_site_motifs = extract_every_rzbs(downstream_seq, template.stem_lens['S3'], template.rzbs_len, template.rzbs_extension_lens['LIVE'])
    num_cut_site_motifs = len(rzbs_cut_site_motifs)
    if num_cut_site_motifs == 0:
        errors['cut_site_motifs'].append("No potential cut sites detected. Unable to process request.")
        return None, errors
    
    data = {'rzbs_positions':rzbs_positions, 'rzbs_cut_site_motifs':rzbs_cut_site_motifs, 'num_cut_site_motifs':num_cut_site_motifs}
    return data, None   

def validate_selected_cut_site_motifs(rzbs_positions, rzbs_cut_site_motifs, selected_cut_site_motifs):
    errors = {'cut_site_motifs':[]}

    if not selected_cut_site_motifs:
        errors['cut_site_motifs'].append('Please select at least one cut site motif')
        return None, errors

    allowed_rzbs_positions = []
    for (position, cut_site_motif) in zip(rzbs_positions, rzbs_cut_site_motifs):
        if cut_site_motif in selected_cut_site_motifs:
            allowed_rzbs_positions.append(position)
    if len(allowed_rzbs_positions) == 0:
        errors['cut_site_motifs'].append('The sum of the counts of the selected cut site motifs must exceed 0')
        return None, errors
        
    data = {'allowed_rzbs_positions':allowed_rzbs_positions, 'selected_cut_site_motifs':selected_cut_site_motifs}
    return data, None

def validate_environmental_parameters(target_environment, in_vivo_option, in_vitro_option, custom_promoter_text):

    errors = {'target_environment':[], 'in_vivo_option':[], 'in_vitro_options':[]}

    if not target_environment:
        errors['target_environment'].append('Please select a target environment.')
        return None, errors

    if target_environment == 'in_vivo':
        if not in_vivo_option:  # Check if in_vivo_option is empty
            errors['in_vivo_option'].append('Please select an option for In Vivo.')
            return None, errors        
        environment = target_environment + '_' + in_vivo_option

    elif target_environment == 'in_vitro':
        if not in_vitro_option:  # Check if in_vitro_option is empty
            errors['in_vitro_option'].append('Please select an option for In Vitro.')
            return None, errors
            
        environment = target_environment + '_' + in_vitro_option
        if in_vitro_option == 'custom_promoter':
            valid_chars = set("AUGTCaugtc-")
            custom_promoter_seq = custom_promoter_text.replace(" ", "")
            if len(custom_promoter_seq) == 0:
                errors['in_vitro_option'].append("Please provide a sequence for your custom promoter")
                return None, errors
            if not set(custom_promoter_seq).issubset(valid_chars):
                errors['in_vitro_option'].append("Custom promoter contains invalid characters. Allowed characters are A, U, G, T, C")
                return None, errors

    data = {'environment':environment}
    return data, None

def validate_temperature(temperature):
    
    errors = {'temperature':[]}

    if not temperature:
        errors['temperature'].append("Please input the environment's temperature")
        return None, errors

    temperature = convert_to_float(temperature) # Returns none if it is not a float
    if temperature == None:
        errors['temperature'].append("Please input a valid temperature")
        return None, errors
    elif temperature < config.min_temp or temperature > config.max_temp:
        errors['temperature'].append("Temperature must be between " + str(config.min_temp) + ' and ' + str(config.max_temp) + ' degrees Celcius.')
        return None, errors
        
    data = {'temperature':temperature}
    return data, None