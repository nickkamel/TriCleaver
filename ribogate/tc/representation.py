import numpy as np

from ..common.utility import arg_sort, inv_list, prefix_sum
from ..common.representation import determine_valid_nucleotide_assignments, build_bcc_to_scc_map, calculate_mutation_probabilities, enforce_bcc_constraints

def update_representation(template, rzbs_position, old_sccs, is_post_processed):
    ''' 
    Updates an individual's segment connected components (SCCs), bcc_to_scc map, and mutation probabilities
    This should be called when the RzBS is mutated or when we are doing post processing with a longer RzBS extension
    '''

    '''
    Keep track of the individual's previous generation bccs (one bcc list per scc object)
    We need these for later because the sccs returned by template.build_sccs have empty bccs
    '''
    old_bccs_list = [scc.bccs for scc in old_sccs] 
    
    # Build the list of segment connected components for the current ribozyme binding site (RzBS)
    sccs = template.build_sccs(rzbs_position, is_post_processed) 

    # Determine the valid nucleotides that the bccs in every scc can assume
    for scc in sccs:
        scc.valid_nucs = determine_valid_nucleotide_assignments(scc.constraints, scc.mismatches, scc.length, scc.size)           

    '''
    Build the bcc_to_scc map that will be used during mutation
    This map depends only on the length of SCCs, which are constant throughout the EA
    However, during the post-processing stage different individuals can have different upstream and downstream rzbs lenths 
    Since the RzBS is entirely contained in the sequence downstream of the repeats, upstream RzBS extensions near the 5' end of the downstream sequence will be truncated
    Same is true for downstream RzBS extensions near the 3' end of the donwstream sequence    
    '''    
    scc_lengths = [scc.length for scc in sccs]
    bcc_to_scc = build_bcc_to_scc_map(scc_lengths)
    
    # Calculate the mutation probabilities
    mutation_probabilities = calculate_mutation_probabilities(bcc_to_scc, sccs)

    '''
    For all non-RzBS segments, the constraints are constant throughout the EA
    However, the constraints of the RzBS segments depend on the current RzBS which varies during the EA and between individuals
    Therefore, we may need to modify certain bccs to ensure that they respect any new RzBS constraints that may occur
    '''
    for i, scc in enumerate(sccs):  
        scc.bccs = enforce_bcc_constraints(old_bccs_list[i], scc.length, scc.valid_nucs)

    return sccs, bcc_to_scc, mutation_probabilities


def create_representation(template, rzbs_position):
    ''' 
    Creates an individual's iniital segment connected components (SCCs), bcc_to_scc map, and mutation probabilities
    '''

    sccs = template.build_sccs(rzbs_position, False) # This clears all previous bccs if there were any

    for scc in sccs:
        scc.valid_nucs = determine_valid_nucleotide_assignments(scc.constraints, scc.mismatches, scc.length, scc.size)           

    # Build mutation map. This should be the same throughout the EA since SCCs have fixed length. But it can change during post-processing.
    scc_lengths = [scc.length for scc in sccs]
    bcc_to_scc = build_bcc_to_scc_map(scc_lengths)

    mutation_probabilities = calculate_mutation_probabilities(bcc_to_scc, sccs)

    return sccs, bcc_to_scc, mutation_probabilities
