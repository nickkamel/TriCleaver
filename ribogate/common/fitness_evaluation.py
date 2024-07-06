import numpy as np

from .utility import avg
import config

def calculate_stem_score(bppm, base_ids_a, base_ids_b):
    base_ids = list(zip(base_ids_a, base_ids_b))
    stem_disruption_score = 0
    for base_id, partner_base_id in base_ids:
        paired_prob = bppm[(base_id, partner_base_id)]
        stem_disruption_score += 1 - paired_prob
    stem_disruption_score /= len(base_ids)
    return stem_disruption_score    

# This function assigns a score from 0 (no disruption) to 1 (full disruption) to each motif
# In this case for TC, the motifs are C1, C2, C3, stem1, stem2, stem3. In TS C1, C2, C3 are treated as a single core motif
def evaluate_rz_activity(bppm, name_to_bases, ribozyme_type):
    
    if ribozyme_type not in ['minimal_hammerhead', 'extended_hammerhead']:
        raise NotImplementedError(f"There is no fitness function for ribozyme type ({ribozyme_type})")    
    is_extended_hammerhead = ribozyme_type == 'extended_hammerhead'

    # Core
    if not config.should_treat_core_as_single_motif:
        #Old version used to generated sRzs. Treats each core segment as its own motif
        core_scores = list()
        for seg_name in ['C1', 'C2', 'C3']:
            seg_score = 0
            base_ids = name_to_bases[seg_name]
            for base_id in base_ids:
                paired_prob = 1 - bppm[(base_id, None)]
                seg_score += paired_prob
            seg_score /= len(base_ids)
            core_scores.append(seg_score)
            
    else:
        # Alternative version used by TS. This will now be used for all future TC designs as well
        # Calculates a single motif score for the entire core. I don't expect the different core scoring methods to noticeably affect EA convergence
        base_ids = name_to_bases['C1'] + name_to_bases['C2'] + name_to_bases['C3']
        core_score = 0
        for base_id in base_ids:
            paired_prob = 1 - bppm[(base_id, None)]
            core_score += paired_prob
        core_score /= len(base_ids)
        core_scores = [core_score]

    # Pseudoknot for extended hammerhead ribozyme. We can't directly predict the pseudoknot, but it is more likely to form if the two knots segments are not paired anywhere else on the ribozyme
    if is_extended_hammerhead:
        knot_score = 0
        base_ids = name_to_bases['S1KNOT'] + name_to_bases['S2KNOT']
        for base_id in base_ids:
            paired_prob = 1 - bppm[(base_id, None)]
            knot_score += paired_prob
        knot_score /= len(base_ids)


    stem_scores = list()
    for i in range(1,4):
        # For the extended hammerhead, the first stem contains a bulge loop which 'splits' the stem into two parts. We term those two parts S0 and S1 (this is a non-official naming convention)
        if is_extended_hammerhead and i == 1:
            base_ids_a = name_to_bases['S0A'] + name_to_bases['S1A']
            base_ids_b = (name_to_bases['S1B'] + name_to_bases['S0B'])[::-1]
        else:
            stem_name = 'S' + str(i)
            base_ids_a = name_to_bases[stem_name + 'A']
            base_ids_b = name_to_bases[stem_name + 'B'][::-1]
        stem_score = calculate_stem_score(bppm, base_ids_a, base_ids_b)
        stem_scores.append(stem_score) 

    if is_extended_hammerhead:
        rz_motif_scores = core_scores + stem_scores + [knot_score]
    else:
        rz_motif_scores = core_scores + stem_scores

    on_score  = 1 -avg(rz_motif_scores) # This reaches max value when the core and stem motifs experience no disruption
    off_score = max(stem_scores) # This reaches max value when at least one stem is completely disrupted

    return on_score, off_score

def calculate_on_off_viability_scores(bppms, truth_vector, name_to_bases, ribozyme_type):
    on_vector = list()
    off_vector = list()

    num_states = len(bppms)
    for i in range(num_states):       
        state_on_score, state_off_score = evaluate_rz_activity(bppms[i], name_to_bases, ribozyme_type) 

        if truth_vector[i] == 1:
            on_vector.append( state_on_score )
        else:
            off_vector.append( state_off_score )

    if len(on_vector) > 0: 
        avg_on = avg(on_vector)
        min_on = min(on_vector)
    # In TS, if the user specifies a functio that is always OFF, the on_vector will be empty.
    # To handle this, we give max ON scores (1)   
    else: 
        avg_on = 1
        min_on = 1
    on_score = avg([min_on, avg_on])

    if len(off_vector) > 0:
        avg_off = avg(off_vector)    
        min_off = min(off_vector)
    # In TS, if the user specifies a functio that is always ON, the off_vector will be empty.
    # To handle this, we give max OFF scores (1)  
    else:
        avg_off = 1
        min_off = 1
    off_score = avg([min_off, avg_off])
    viability_score = avg([min_on, min_off]) # Has a value between 0 and 1

    return on_score, off_score, viability_score

def coarse_grain_structures(bppms, base_to_seg, extension_region_seg_ids, truth_vector, are_segments_low_level):
    '''

    Coarse grains multiple base-pairing probability matrices (BPPMs) into their corresponding segment-pair magntiude matrices (SPMMs)
    The switchable ribozyme can be viewed as two components: the ribozyme (Rz) and the extension region (ER).
    We can subdivide the extension region into more fine-grained segments while keeping the ribozyme as a single abstract segment (Rz). This is called a high-level (H) segment partitioning.
    For example, in TS the H segments may be: Rz, OBS1, OBS2. In TC they could be Rz, L0, SENSOR, L1.
    Alternatively, we can subdivide both the extension region and the ribozyme into fine-grained segments. This is called a low-level (L) segment partitioning.
    For example, in TS the L segments could be Stem1A, C1, Stem2A, OBS1, OBS2, Stem2B, C2, etc.
    An SPPM encodes the relative weight of any two L or H segments being paired together.

    [TO DO: In Summer 2024 I added pre-RzBS and post-RZBS segments for TC. These are neither extension region segments nor part of the ribozyme. 
            This was not a scenario I originally accounted for. With the way the code is currently written, pre-RzBS and post-RzBS are treated as ribozyme segments.
            I don't think this breaks anything (especially for L segments), but we might want to keep an eye on this just to be sure.]
    
    Arguments:
        bppms (list of dicts): One base-pairing probability matrix per state. bppm[i,j] stores the probability of nucleotides i and j being paired. bppm[i, None]  stores the probability of nucleotide i being unpaired
        base_to_seg (dict). Map from integer global nucleotide position to pair (position of the segment within the strand, position of base within that segment) 
        extension_region_seg_ids (list of ints): The positions of the extension regon segments in the strand
        truth_vector (boolean list): The target ribozyme condition (0 = inactive, 1 = active) in each state
        are_segments_low_level (bool): If true coarse-grains using L-segments; if false uses H-segments. 

    Returns:
        spmms (mun_states X num_segs X num_segs numpy float array): One segment pair magnitude matrix for each state
    '''

    if are_segments_low_level:
        # Every segment is considered
        num_segs = max( [v[0] for v in base_to_seg.values()]) + 1 
    else:
        # Only consider the extension region segs + an abstract segment consisting of the entire ribozyme
        num_segs = len(extension_region_seg_ids) + 1
        er_start_seg_id = min(extension_region_seg_ids)

    num_states = len(bppms)
    spmms = np.zeros(( num_states, num_segs, num_segs ))
    for state_id in range(num_states):
        for nucleotide_pair_key, pairing_prob in bppms[state_id].items(): # [Old comment not sure if still relevant]: dict.items() is super slow if bppm is full
            if None not in nucleotide_pair_key: # Ignore unpaired bases
                
                # Find the segment of each nucleotide in the pair
                seg1_id = base_to_seg[nucleotide_pair_key[0]][0]
                seg2_id = base_to_seg[nucleotide_pair_key[1]][0]
                
                if are_segments_low_level:
                    if truth_vector[state_id] == 1:
                        # When the target output is 1, only consider interactions between the extension region segments.
                        if seg1_id in extension_region_seg_ids and seg2_id in extension_region_seg_ids:
                            spmms[state_id][seg1_id][seg2_id] += pairing_prob
                    else: # In states that output a 0, consider all possible segment interactions
                        spmms[state_id][seg1_id][seg2_id] += pairing_prob

                else:
                    '''
                    If we coarse-grain using high-level (H) segments, we consider two types of interactions:
                    1) An extension region (ER) segment and another ER segment. These are considered in each state.
                    2) An ER segment and the coarse-grained ribozyme segment. These are only considered in states that output a 0.
                    
                    Instead of keeping track of the specific ribozymes segments (e.g. Stem1A, CoreA, Stem2A, etc.), we treat the ribozyme as a single Rz segment at position 0
                   
                     Given N ER segments, we consider N ER-Rz and N Rz-ER interactions (these are symmetric).
                     .... In other words we treat the ribozyme 
                    Rz-Rz interactions are never considered
                    '''

                    # If at least one of the segs is an extension region (ER) segment
                    if seg1_id in extension_region_seg_ids or seg2_id in extension_region_seg_ids:
                        # Always consider ER-ER interactions regardless of target output
                        if seg1_id in extension_region_seg_ids and seg2_id in extension_region_seg_ids:
                            spmms[state_id][seg1_id - er_start_seg_id + 1][seg2_id - er_start_seg_id + 1] += pairing_prob
                            
                        '''
                        Only consider ER-Rz and Rz-ER interactions if the target output is 0.
                        When accumulating the probabilities, we don't keep track of the specific Rz segment that the ER paired with. 
                        We just treat it as interacting with a coarse ribozyme segment at position 0
                        '''
                        if truth_vector[state_id] == 0:
                            if seg1_id in extension_region_seg_ids and seg2_id not in extension_region_seg_ids: # ER-Rz
                                spmms[state_id][seg1_id - er_start_seg_id + 1][0] += pairing_prob # Using [0] to store Rz seg vals. This is different from paper which uses [num_segs - 1]
                            if seg1_id not in extension_region_seg_ids and seg2_id in extension_region_seg_ids: # Rz-ER
                                spmms[state_id][0][seg2_id - er_start_seg_id + 1] += pairing_prob

    return spmms # This is symmetric since in Vienna, we make the bppm dict symmetric