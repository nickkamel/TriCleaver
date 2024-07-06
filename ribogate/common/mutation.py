import numpy as np

def mutate_bcc(bcc_to_scc, mutation_prob, sccs):
    '''
    Mutates a base connected component (BCC), i.e. mutates a single nucleotide or a pair of complementary nucleotides
        
    Parameters:
        bcc_to_scc (dict): A map from the index of a bcc to its position within its parent SCC (segment connected component)
        mutation_prob (list): The probability distrubtion of a given bcc being selected for mutation
        sccs (list of SCC objects): The individual's SCC before mutation
            
    Returns
        sccs (list of SCC objects): The individual's SCCs after mutation            
    '''
        
    num_bccs = len(bcc_to_scc)
    # Randomly select a (global) index of a bcc and use the bcc_to_scc map to get the scc containing that bcc and the bcc's relative position within that scc
    bcc_id_global = np.random.choice( range( num_bccs),  1, p=mutation_prob)[0] 
    (scc_id, bcc_id_scc) = bcc_to_scc[bcc_id_global] 
    scc = sccs[scc_id]
    bcc = scc.bccs[bcc_id_scc]
        
    # Change the bcc to a different valid bcc
    options = [i for i in scc.valid_nucs[bcc_id_scc] if i != bcc]                              
    scc.bccs[bcc_id_scc] = np.random.choice(options)      
        
    return sccs
