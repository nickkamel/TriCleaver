import numpy as np
import copy

from ..common.folding_constraints import single_nucleotide_constraint_lut, paired_nucleotide_constraint_lut
from ..common.utility import prefix_sum, arg_sort

class SwitchableRibozyme():
    def __init__(self):
        '''
        Represents individuals evolved by TruthSeqEr and TriCleaver
        Contains no methods: it just being used as a struct
        '''

        '''
        SCCs: Segment Connected Components objects that enforce sequence constraints and structural constraints for unpaired segments or between segment pairs
        segments: Segment objects that are generated from the SCCs
        strands: The segments are finally merged into complete strands. Each strand is a string of characters from {A, U, G, C}
        '''
        self.sccs = list()
        self.segments = list()
        self.strands = list()

        '''
        Map used for generated folding constraints and evaluating ribozyme activity
        input: name of a segment
        output: list of base (global) positions contained in that segment
        Here 'global' means that the nucleotides in a set of multiple strands are numbered consecutively and don't reset        
        Example: Strand 0 = 5' AAAAA 3' and strand 1 = GGGGG. The first nucleotide of strand 1 has global position 5 since it is preceded by 5 nucleotides on strand 0.
        '''
        self.name_to_bases = dict() 
        
        '''
        Map used to coarse-grain the secondary structures adopted by the switchable ribozyme
        input: global id of a given base
        output: (global id of a the segment containing that base, position of that base within that segment)
        Here 'global' means that the segments in a set of multiple strands are numbered consecutively
        '''
        self.base_to_seg = dict()  

        '''
        Map used for mutation
        Each segment connected component (SCC) can be broken down into a list of base connected components (BCC)
        These bccs can be merged into a single list and assigned a global id
        input: global id of a base connected component (BCC)
        output: (id of the SCC containing that BCC, position of the BCC within the SCC that contains it)
        '''
        self.bcc_to_scc = dict()
        
        # The probability of a given BCC being selected for mutation
        self.bcc_mutation_probabilities = list()
        
        # The index of the first nucleotide of the ribozyme binding site w.r.t. sequence downstream of the repeats region
        self.rzbs_position = None
        
        # A list of (sequence, folding constraints) to allow the folding software to simulate the switchable ribozyme binding to different combinations of inputs or repat lengths.
        self.tasks = list()

        '''
        The evaluated performance of the switchable ribozyme across a variety of metrics.
        Each value can be potentially used for parent/survivor selection, but not all of them need to be.
        '''
        self.potential_fitnesses = dict()

class SegmentConnectedComponent():
    def __init__(self, name, strands, positions, constraints, mismatches=[]):
        self.name = name
        self.size = len(name) # 1 for unpaired, 2 for paired
        self.strands = strands
        self.positions = positions
        if self.size == 2:
            self.constraints = list(zip(constraints[0], constraints[1][::-1]))
        elif self.size == 1:
            self.constraints = list(constraints[0])

        self.length = len(self.constraints)
        self.mismatches = mismatches
        self.bccs = list()
             
class Segment():
    def __init__(self, name, length, strand, pos, seq):
        self.name = name
        self.length = length
        self.strand = strand
        self.pos = pos
        self.seq = seq

################

def determine_segment_positions(segment_names):
    '''
    Determines the positions of an ordered list of segments with respect to the 5' end of the strand
        
    Parameters:
        segment_names (list of strings)
            
    Returns:
        segment_positions (dict). Maps segment names to their (0-indexed) position w.rt. to 5' end of strand
        
    '''
    
    segment_positions = dict()
    position = 0
    for name in segment_names:
        segment_positions[name] = position
        position += 1
            
    return segment_positions

def initialize_bccs(length, valid_nucs):
    '''
    For a segment connect component (SCC) of a certain length and with certain constraints, randomly samples a list of valid base connected components (BCCs)
    '''
    
    bccs = list()
    for i in range(length):
        bccs.append( np.random.choice(valid_nucs[i]))
    return bccs


def enforce_bcc_constraints(old_bccs, new_length, valid_nucs):
    old_length = len(old_bccs)           
    bccs = copy.deepcopy(old_bccs)
            
    # If a bcc becomes invalid due to a change in constraints, assign that bcc a new random valid value    
    for i in range(old_length):
        if bccs[i] not in valid_nucs[i]:
            bccs[i] = np.random.choice(valid_nucs[i])
       
    # Handle length changes
    # Add random bccs if the new length is greater than the old one        
    if new_length > old_length:
        for i in range(old_length, new_length):
            bccs.append( np.random.choice(valid_nucs[i]) ) 

    # Truncate the bccs if the new length is smaller than the hold one
    elif new_length < old_length:
        bccs = bccs[0:new_length]
                
    return bccs
                   

def determine_valid_nucleotide_assignments(constraints, mismatches, length, size):
    '''
    Based on the SCC's constraints, generates a set of valid nucleotide assignments

    Parameters:
        constraints (list of pairs / chars): If the SCC is of size 1 (unpaired), each element in the constraints list is a single character representing the constraint of an unpaired nucleotide. If the SCC is of size 2 (paired), each element is a pair representing the pairs of constraints of a pair of nucleotides
        mismatches (list of ints): A list of positions where we want to the nucleotides to be non-complementary. Only applicable to SCCs of size 2
        length (int): Number of elements (bases / base-pairs) in the SCC
        size (int): Whether the SCC represented an isolated segment (size 1) or a pair of reverse complementary segments (size 2)

    Returns:
        valid_nucs (list of list of strings): List of allowed nucleotide / nucleotide pair characters for each nucleotide / nucleotide pair in the constraint 
    '''

    valid_nucs = list()
    for i in range( length ):
        if size == 1: #unpaired                
            cur_valid = single_nucleotide_constraint_lut( constraints[i] )
        else:
            is_mismatch = i in mismatches
            cur_valid = paired_nucleotide_constraint_lut( constraints[i][0], constraints[i][1], is_mismatch )
        valid_nucs.append( cur_valid )
    return valid_nucs

def build_bcc_to_scc_map(lengths):
    '''
    Builds a positional map from nucleotide connected components to segment connected components to be used by the mutation operator.

    Parameters:
        lengths (list of ints): The number of nucleotide connected components in each segment connected component

    Returns:
        bcc_to_scc (dict): map from a 0-indexed scalar integer representing the global id of a bcc to a pair of integers representing the scc that contains it and its local id within that scc
    '''
    
    num_sccs = len(lengths)

    bcc_to_scc = dict()
    p_sum = [0] + prefix_sum(lengths)

    for scc_id in range(num_sccs):
        # local_id is the position of the ncc within the scc
        # global_id is the position of the ncc when the sccs are treated as a single contiguous element
        for local_id, global_id in enumerate( range(p_sum[scc_id], p_sum[scc_id+1]) ):
            bcc_to_scc[global_id] = (scc_id, local_id) 
    return bcc_to_scc

def calculate_mutation_probabilities(bcc_to_scc, sccs):
    '''
    Calculates the probability that a given base connected component (BCC) will be selected for mutation
    '''

    num_bccs = len(bcc_to_scc)
    weights = list()
    for i in range(num_bccs):
        
        # For a given BCC, the determine the SCC it is located in and its position with that SCC
        (scc_id, bcc_id_loc) = bcc_to_scc[i]
        scc = sccs[scc_id]
        
        # Assign that BCC a weight depending on the number of valid nucleotide assignment is had
        cur_valid = scc.valid_nucs[bcc_id_loc]
        if scc.size == 1: # unpaired
            weight = len(cur_valid) - 1  # -1 because we care about the #of valid bccs that a bcc can mutate into 
        else:
            weight = ( len(cur_valid) - 1 ) / float(2) #/2 because each base pair contains 2 bases
        weights.append(weight)
        
    # Normalize the weights so that they sum to 1
    denom = sum(weights)
    bcc_mutation_probabilities = [w/float(denom) for w in weights]

    return bcc_mutation_probabilities
