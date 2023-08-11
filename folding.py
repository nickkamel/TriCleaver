from utility import prefix_sum

def allowed_lut(constraint):
    #form allowed LUT
    #Obviously forming the LUT each time the function is called is ridiculous
    #But if I form it once and store it in EA, inds need to acess it, and gens need to acess it from inds
    allowed_lut = dict()
    allowed_lut[('N', 'N', False)] = ['GC', 'CG', 'GU', 'UG', 'AU', 'UA']
    allowed_lut[('N', 'N', True )] = ['GA', 'AG', 'GG', 'UU', 'CC', 'AA', 'AC', 'CA', 'UC', 'CU']
    allowed_lut[('G', 'N', False)] = [bp for bp in allowed_lut[('N', 'N', 0)] if bp[0] == 'G'] 
    allowed_lut[('H', 'N', False)] = [bp for bp in allowed_lut[('N', 'N', 0)] if bp[0] != 'G']
    allowed_lut[('N', 'H', False)] = [bp for bp in allowed_lut[('N', 'N', 0)] if bp[1] != 'G'] 
    allowed_lut[('G', 'C', False)] = [('GC')] 
    allowed_lut[('G', 'U', False)] = [('GU')]
    allowed_lut[('C', 'G', False)] = [('CG')] 
    allowed_lut[('U', 'G', False)] = [('UG')] 
    allowed_lut[('U', 'A', False)] = [('UA')] 
    allowed_lut[('A', 'U', False)] = [('AU')] 
    allowed_lut[('N', 'C', False)] = [('GC')] 
    allowed_lut[('N', 'U', False)] = [('GU'), ('AU')] 
    allowed_lut[('N', 'A', False)] = [('UA')] 
    allowed_lut[('N', 'G', False)] = [('UG'), ('CG')] 
    allowed_lut[('N', 'G', True )] = [('AG'), ('GG')] 
    allowed_lut[('N', 'C', True )] = [('AC'), ('CC'), ('UC')] 
    allowed_lut[('N', 'A', True )] = [('AA'), ('CA'), ('GA')] 
    allowed_lut['A']               = ['A']
    allowed_lut['U']               = ['U']
    allowed_lut['G']               = ['G']
    allowed_lut['C']               = ['C']
    allowed_lut['N']               = ['A', 'U', 'G', 'C']
    allowed_lut['H']               = ['A', 'U', 'C']

    return allowed_lut[constraint]

def build_map(lengths, bool1_idx): #whether the global coordinate system is 1-indexed or 0-indexed
    map_f  = dict()
    map_f_h = list() #Hierarchical
    map_r  = dict()
    p_sum  = prefix_sum(lengths)
    p_sum0 = [0] + p_sum
    if bool1_idx == 1:
        p_sum0 = [p + 1 for p in p_sum0]

    for i in range( len(p_sum) ): #for each segment
        seg_ids_global = list() #FH
        for local_id, global_id in enumerate( range(p_sum0[i], p_sum0[i + 1]) ):
            map_f[(i, local_id)] = global_id #Foward map #ith segment, "local" ith base within segment
            map_r[global_id]     = (i, local_id)  #Reverse map
            seg_ids_global.append(global_id) #FH 
        map_f_h.append(seg_ids_global) #FH
    return (map_f_h, map_f, map_r)


def generate_fold_constraint(seqLen, constrainedBaseIds):
    constraints = ''
    for b in range(1, seqLen + 1): #We add the +1 here because the constrainedBaseIds are 1-indexed
        if b in constrainedBaseIds:
            constraints += 'x'
        else:
            constraints += '.'
    return constraints





