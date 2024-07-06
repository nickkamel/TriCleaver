import copy

from ..common.representation import SegmentConnectedComponent, determine_segment_positions
import config

class RibozymeTemplate:  
    def __init__(self, ribozyme_type, template_params, downstream_sequence, repeats_sequence):
        '''
        Each individual is represented by a state dependency graph (SDG).
        The SDG is a set of segmment connected components (SCC).
        Each SCC consists of one or two segments and stores the positions of these segments in the sRz / substrate strand as well as sequence constraints. 
        The segments in the sRz are fixed throughout the EA.
        However, those in the ribozyme binding site (RzBS) of the substrate will vary throughout the EA since the RzBS is can change
        Therefore, the full SDG can only be completed once the RzBS is known
        Here, we specify the part of the SCCs that is constant throughout the EA
        
        Parameters:
            ribozyme_type (string): Name of the type of ribozyme we are using
            downstream_sequence (string): The sequence of the substrate that is downstream of the repeats
            repeats_sequence (string): The sequence of repeats        
        '''
        
        self.ribozyme_type = ribozyme_type
        self.downstream_sequence = downstream_sequence
        self.repeats_sequence = repeats_sequence
        
        # Specify segment names and lengths. Specify which segments are connected together in segment connected components (SCCs)
        self.sensor_len = len(repeats_sequence)  
        self.linker_len = template_params['LINKER_LEN']
        self.rzbs_extension_lens = template_params['RZBS_EXTENSION_LENS']
        if self.ribozyme_type == 'minimal_hammerhead':
            self.scc_names = [('S1A', 'S1B'), ('C1',), ('S2A', 'S2B'), ('L0',), ('SENSOR','REPEATS'), ('L1',), ('C2',), ('S3A', 'S3B'), ('UPSTREAM_RZBS',), ('C3',), ( 'DOWNSTREAM_RZBS',)] # Each name is a tuple, not a list so it can serve as dict key later
            self.srz_segment_names = ['S1A', 'C1', 'S2A', 'L0', 'SENSOR', 'L1', 'S2B', 'C2', 'S3A'] # srz: selective ribozyme
            self.rzbs_segment_names = ['UPSTREAM_RZBS', 'S3B', 'C3', 'S1B', 'DOWNSTREAM_RZBS'] # rzbs: ribozyme binding site 
            
            self.stem_lens = template_params['STEM_LENS_MINIMAL']
            self.rzbs_len = self.stem_lens['S1'] + 1 + self.stem_lens['S3'] # +1 is due to the single nucleotide in the third part of the core

        elif self.ribozyme_type == 'extended_hammerhead':
            '''
            Stem 1 of the extended hammarhead contains a bulge loop. 
            To handle this, I divide stem 1 into two parts flanking the loop. 
            For lack of better terminology, I refer to these two halves as stem 0 and stem 1.
            '''
            self.scc_names = [('S0A', 'S0B'), ('S1A', 'S1B'), ('S1KNOT',), ('C1',), ('S2A', 'S2B'), ('L0',), ('S2KNOT',), ('SENSOR','REPEATS'), ('L1',), ('C2',), ('S3A', 'S3B'), ('UPSTREAM_RZBS',), ('C3',), ( 'DOWNSTREAM_RZBS',)]
            self.srz_segment_names = ['S0A', 'S1KNOT', 'S1A', 'C1', 'S2A', 'S2KNOT', 'L0', 'SENSOR', 'L1', 'S2B', 'C2', 'S3A'] 
            self.rzbs_segment_names = ['UPSTREAM_RZBS', 'S3B', 'C3', 'S1B', 'S0B', 'DOWNSTREAM_RZBS']         
           
            self.stem_lens = template_params['STEM_LENS_EXTENDED']
            self.rzbs_len = self.stem_lens['S0'] + self.stem_lens['S1'] + 1 + self.stem_lens['S3']
       
        self.specify_segment_positions()
        self.specify_segment_constraints()

        # Prepare the sensor mismtaches
        self.num_contiguous_mismatches = template_params['NUM_CONTIGUOUS_MISMATCHES']
        self.mismatch_spacing = template_params['MISMATCH_SPACING']
        self.scc_mismatches = {key: [] for key in self.scc_names if len(key) == 2} # Set default for mismatches (which can only potentially be defined for SCCs of size 2)
        self.scc_mismatches[('SENSOR','REPEATS')] = self.place_mismatches(self.num_contiguous_mismatches, self.mismatch_spacing, self.sensor_len)
        
        # Specify max # of strides. Each stride is a different way in which we simulate the wild-type repeats binding to the sensor
        self.max_num_strides = template_params['MAX_NUM_STRIDES']
        
    def specify_segment_positions(self):
        self.segment_strands = dict()
        self.segment_positions = dict()
        
        srz_segment_positions = determine_segment_positions(self.srz_segment_names)
        for name in self.srz_segment_names:
            self.segment_strands[name] = 0
            self.segment_positions[name] = srz_segment_positions[name]

        rzbs_segment_positions = determine_segment_positions(self.rzbs_segment_names)
        for name in self.rzbs_segment_names:
            self.segment_strands[name] = 1
            self.segment_positions[name] = rzbs_segment_positions[name]

        # The repeats are on a virtual 3rd strand that is not co-folded. It is used to constrain the sensor sequence. Its effect is taken into account using a folding constraint string.
        self.segment_strands['REPEATS'] = 2
        self.segment_positions['REPEATS'] = 0

        self.extension_region_seg_ids = list()
        for i, name in enumerate(self.srz_segment_names):
            if name in ['L0', 'SENSOR', 'L1']:
                self.extension_region_seg_ids.append(i)   
                
        self.num_strands = 3

    def specify_segment_constraints(self):
        self.segment_constraints_const = dict() 
        
        if self.ribozyme_type == 'extended_hammerhead':
            self.segment_constraints_const['S0A'] = 'N'*self.stem_lens['S0']
            self.segment_constraints_const['S1KNOT'] = 'AAU'
            self.segment_constraints_const['S2KNOT'] = 'UGAAAU'

        self.segment_constraints_const['S1A'] = 'N'*self.stem_lens['S1']     
        self.segment_constraints_const['C1'] = 'CUGAUGA'
        self.segment_constraints_const['S2A'] = 'G' + 'N'*(self.stem_lens['S2']-1)
        self.segment_constraints_const['L0'] = 'N'*self.linker_len
        self.segment_constraints_const['SENSOR'] = 'N'*self.sensor_len
        self.segment_constraints_const['L1'] = 'N'*self.linker_len
        self.segment_constraints_const['S2B'] = 'N'*(self.stem_lens['S2']-1) + 'C'
        self.segment_constraints_const['C2'] = 'GAA'
        self.segment_constraints_const['S3A'] = 'A' + 'N'*(self.stem_lens['S3']-1)

        # All the constraints below depend on the RzBS which changes during the EA run
        for name in self.rzbs_segment_names:
            self.segment_constraints_const[name] = None

        self.segment_constraints_const['REPEATS'] = self.repeats_sequence
            
    def place_mismatches(self, num_contiguous_mismatches, spacing, segment_len):
        '''
        Generates a list of contiguous mismatches groups evenly spaced out along a segment
        These mismatches give ribozyme some flexiblity when bound to transcript. 
        They also avoid extensive dsRNA that could potentially cause PKR or Dicer problems

        Parameters:
            num_contiguous_mismatches (int): the size of a contiguous group of mismatches
            spacing (int): number of nucleotides between the mismatches. If 0, mismatches will be consecutive
            segment_len (int): the length of the segment for which we are placing the mismatches

        Returns:
            mismatch_idxs (list of ints): The 0-indexed locations of the mismatches on the segment
        '''

        mismatch_idxs = list()
        group_start_idx = spacing + 1 # The +1 is just to be backwards compatible with my old approach
        while group_start_idx < segment_len:
            group_end_idx = min(group_start_idx + num_contiguous_mismatches, segment_len) # Truncates the last group to avoid it going out of bounds
            for i in range(group_start_idx, group_end_idx):
                mismatch_idxs.append(i)
            group_start_idx += num_contiguous_mismatches + spacing

        return mismatch_idxs
        
    def extract_rzbs_segments(self, rzbs_position, is_post_processed):
        '''
        Extracts the sequences of the RzBS segments for a given RzBS position
        
        Parameters:
            rzbs_position (int): The index of the first nucleotide of the ribozyme binding site w.r.t to the start of the downstream sequence
            is_post_processed (bool): If True indicates that we are in post-processing statge of EA; if False indicates that we are live in the EA loop
        Returns
            rzbs_segment_sequences (dict of strings): The sequence of each RzBS segment        
        '''

        # If we are in post-processing stage, use the POST RzBS extension length (which may be longer than the LIVE length)
        if is_post_processed:
            rzbs_extension_len = self.rzbs_extension_lens['POST'] 
        else:                        
            rzbs_extension_len = self.rzbs_extension_lens['LIVE']

        rzbs_segment_sequences = dict()
        rzbs_segment_sequences['UPSTREAM_RZBS'] = self.downstream_sequence[max(0, rzbs_position-rzbs_extension_len):rzbs_position] # max(0, ...) ensures that the longer post-processing extension doesn't exceed downstream seq bounds
        rzbs = self.downstream_sequence[rzbs_position:rzbs_position + self.rzbs_len]
        
        s3b_start = 0
        c3_start = self.stem_lens['S3']
        s1b_start = c3_start + 1
        rzbs_segment_sequences['S3B'] = rzbs[0:c3_start]
        rzbs_segment_sequences['C3'] = rzbs[c3_start:s1b_start]
        if self.ribozyme_type == 'minimal_hammerhead':
            rzbs_segment_sequences['S1B'] = rzbs[s1b_start:self.rzbs_len]
        elif self.ribozyme_type == 'extended_hammerhead':
            s0b_start = s1b_start + self.stem_lens['S1']
            rzbs_segment_sequences['S1B'] = rzbs[s1b_start:s0b_start]
            rzbs_segment_sequences['S0B'] = rzbs[s0b_start:self.rzbs_len]
        
        rzbs_segment_sequences['DOWNSTREAM_RZBS'] = self.downstream_sequence[rzbs_position + self.rzbs_len:min(rzbs_position + self.rzbs_len + rzbs_extension_len, len(self.downstream_sequence) )] # min(...) ensures that the longer post-processing extension doesn't exceed downstream seq bounds

        return rzbs_segment_sequences

    def build_sccs(self, rzbs_position, is_post_processed):
        '''
        BUilds the list of segment connected componnets (SCC) for a switchable ribozyme
        This generated SCCs vary depending on to the current RzBS position
        
        Parameters:
            rzbs_position (int): The index of the first nucleotide of the ribozyme binding site w.r.t to the start of the downstream sequence
            is_post_processed (bool): If True indicates that we are in post-processing statge of EA; if False indicates that we are live in the EA loop
        Returns
            sccs (list of SCC objectss): List of segment connected componnets (SCC)
        '''

        # Extract RzBS segments
        rzbs_segment_sequences = self.extract_rzbs_segments(rzbs_position, is_post_processed)

        # Add the extracted rzbs segments to the segment constraints dict
        segment_constraints = copy.deepcopy(self.segment_constraints_const)
        for (seg_name, seg_constraint) in rzbs_segment_sequences.items():
            segment_constraints[seg_name] = seg_constraint

        # We now have all the information required to instantiate every SCC
        sccs = list()
        for name in self.scc_names:
            if len(name) == 1: # unpaired
                scc = SegmentConnectedComponent(name, [ self.segment_strands[name[0]] ], [ self.segment_positions[name[0]] ], [ segment_constraints[name[0]] ], [])
            elif len(name) == 2: # paired
                scc = SegmentConnectedComponent(name, [ self.segment_strands[name[0]], self.segment_strands[name[1]] ], [ self.segment_positions[name[0]], self.segment_positions[name[1]] ], [ segment_constraints[name[0]], segment_constraints[name[1]] ], self.scc_mismatches[name] )
            sccs.append(scc)

        return sccs

