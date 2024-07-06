def single_nucleotide_constraint_lut(constraint_char):
    '''
    Look up table storing the allowed values a single nucleotide can assume given certain constraints
    '''

    lut = dict()
    lut['A'] = ['A']
    lut['U'] = ['U']
    lut['G'] = ['G']
    lut['C'] = ['C']
    lut['N'] = ['A', 'U', 'G', 'C']
    lut['H'] = ['A', 'U', 'C']

    return lut[constraint_char]

def paired_nucleotide_constraint_lut(constraint_char_0, constraint_char_1, bool_mismatch):
    '''
    Look up table storing the allowed values a pair of nucleotide can assume given certain constraints
    If a mismatch is specified, the values must be non-complementary
    '''

    lut = dict()
    lut[('N', 'N', False)] = ['GC', 'CG', 'GU', 'UG', 'AU', 'UA']
    lut[('N', 'N', True )] = ['GA', 'AG', 'GG', 'UU', 'CC', 'AA', 'AC', 'CA', 'UC', 'CU']
    lut[('G', 'N', False)] = [bp for bp in lut[('N', 'N', 0)] if bp[0] == 'G'] 
    lut[('H', 'N', False)] = [bp for bp in lut[('N', 'N', 0)] if bp[0] != 'G']
    lut[('N', 'H', False)] = [bp for bp in lut[('N', 'N', 0)] if bp[1] != 'G'] 
    lut[('G', 'C', False)] = ['GC'] 
    lut[('G', 'U', False)] = ['GU']
    lut[('C', 'G', False)] = ['CG'] 
    lut[('U', 'G', False)] = ['UG'] 
    lut[('U', 'A', False)] = ['UA'] 
    lut[('A', 'U', False)] = ['AU'] 
    lut[('N', 'C', False)] = ['GC'] 
    lut[('N', 'U', False)] = ['GU', 'AU'] 
    lut[('N', 'A', False)] = ['UA'] 
    lut[('N', 'G', False)] = ['UG', 'CG'] 
    lut[('N', 'G', True )] = ['AG', 'GG'] 
    lut[('N', 'C', True )] = ['AC', 'CC', 'UC'] 
    lut[('N', 'A', True )] = ['AA', 'CA', 'GA'] 

    return lut[(constraint_char_0, constraint_char_1, bool_mismatch)]

def generate_folding_constraint_string(seq_len, constrained_base_positions):
    '''
    Builds a constraint string to be used with Vienna RNA

    Parameters:
        seq_len (int): Length of the sequence to which the constraint applies
        constrained_base_positions (list of ints): 1-indexed positions of the bases to be constrained as unpaired

    Returns:
        constraint (string): String where 'x' denotes corresponding base is forced unpaired and '.' denotes no constraint
    '''

    constraint = ''
    for b in range(1, seq_len + 1):
        if b in constrained_base_positions:
            constraint += 'x'
        else:
            constraint += '.'
    return constraint


