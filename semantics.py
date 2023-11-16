from utility import avg, flatten
from vienna  import fold
from fitness import evaluate_rz_activity
from repr    import Segment
import numpy as np

#The switchable ribozyme is can be viewed as two components: the ribozyme (RZ) and the extension region (ER). #This view is called Opaque-Opaque (OO)
#We can subdivide the extension region into more fine-grained segments while keeping the ribozyme as a single abstract segment. This is called Transparent-Opaque (TO)
#We can sudvide both the extension region and the ribozyme into fine-grianed segments. This is called transparent-transparent (TT)
def coarse_grain_pheno(bppms, map, er_seg_ids, truth_vector, bool_TT):
    if bool_TT:
        num_segs         = max( [v[0] for v in map.values()]) + 1 
    else: #TO
        num_segs         = len(er_seg_ids) + 1
        er_start_seg_id = min(er_seg_ids)


    num_states    = len(bppms)
    t_sppms = np.zeros(( num_states, num_segs, num_segs ))
    for state_id in range(num_states):
        for key, val in bppms[state_id].items(): #!!!! dict.items() is super slow if bppm is full!!!
            if None not in key: #Ignore unpaired bases
                seg1_id = map[key[0]][0]
                seg2_id = map[key[1]][0]
                if bool_TT:
                    if truth_vector[state_id] == 1:
                        #When the target output is 1, only consider interactions between the extension region segments.
                        if seg1_id in er_seg_ids and seg2_id in er_seg_ids:
                            t_sppms[state_id][seg1_id][seg2_id] += val
                    else:
                        t_sppms[state_id][seg1_id][seg2_id] += val

                else:
                    if seg1_id in er_seg_ids or seg2_id in er_seg_ids:
                        #always consider ER-ER interactions regardless of target output
                        if seg1_id in er_seg_ids and seg2_id in er_seg_ids:
                            t_sppms[state_id][seg1_id - er_start_seg_id + 1][seg2_id - er_start_seg_id + 1] += val
                        #only consider ER-Rz and Rz-ER interactions if the target output is 0
                        if truth_vector[state_id] == 0:
                            if seg1_id in er_seg_ids and seg2_id not in er_seg_ids:
                                t_sppms[state_id][seg1_id - er_start_seg_id + 1][0] += val #Using [0] to store Rz seg vals. This is different from paper which uses [num_segs - 1]
                            if seg1_id not in er_seg_ids and seg2_id in er_seg_ids:
                                t_sppms[state_id][0][seg2_id - er_start_seg_id + 1] += val


    return t_sppms #This is symmetric since in Vienna, we make the bppm dict symmetric

def pheno_perf_fits(input_data): #Operates on a batch of selected individual properties.
    tasks, common_data = input_data
    base_to_seg, name_to_bases, er_seg_ids, truth_vector, bool_300, folding_temp, bool_TT, mismatches = common_data
    num_states         = len(truth_vector)
    num_inds_batch     = len(tasks)
    logic_tasks        = flatten([task[0] for task in tasks])
    logic_results      = fold(logic_tasks, bool_300, folding_temp, 2) #2 is for 2 strands so it we use cofold
   
    results            = list()   

    #Prepare for mismatches calculation. Each individual has mismatches at the same location.
    sensor_base_ids      = name_to_bases['SENSOR'] #These are 1-indexed and right inclusive
    sensor_lims          = (sensor_base_ids[0] , sensor_base_ids[-1] + 1) #These are 1-indexed and right exclusive
    #mismatches gives the indices of the sensor mismatches w.r.t. to their local position in sensor whereas mismatches_global gives their position w.r.t. the entire strand
    mismatches_global = [o + sensor_lims[0] for o in mismatches] 

    for ind_idx_batch in range(num_inds_batch):              
        potential_fitnesses = dict()
        bppms               = list()
        mfe_dbs             = list()
        mfe_ads             = list()
        free_energies       = list()
        on_scores           = list()
        off_scores          = list()
        mismatch_score      = 0
      
        for state_id in range(num_states):
            task_id = ind_idx_batch*num_states + state_id
            (mfe_db, bppm, mfe_ad, free_energy) = logic_results[task_id] #Cofold has no ensemble diversity

            on_score, off_score = evaluate_rz_activity(bppm, name_to_bases)            

            if truth_vector[state_id] == 1:
                on_scores.append( on_score )
            else:
                off_scores.append( off_score )
           
            mfe_dbs.append(mfe_db)           
            bppms.append(bppm)
            mfe_ads.append(mfe_ad)
            free_energies.append(free_energy)

            ###Mismatch
            for m in mismatches_global:
                mismatch_score += bppm[(m, None)]
        mismatch_score = mismatch_score / float(len(mismatches_global) )

        potential_fitnesses['avgON']     = avg(on_scores)
        potential_fitnesses['avgOFF']    = avg(off_scores)
        potential_fitnesses['minON']     = min(on_scores)
        potential_fitnesses['minOFF']    = min(off_scores)
        potential_fitnesses['ON']        = avg([ potential_fitnesses['minON'], potential_fitnesses['avgON'] ])
        potential_fitnesses['OFF']       = avg([ potential_fitnesses['minOFF'], potential_fitnesses['avgOFF'] ])
        potential_fitnesses['switch']    = potential_fitnesses['ON'] + potential_fitnesses['OFF']
        potential_fitnesses['switchMin'] = sum([ potential_fitnesses['minON'], potential_fitnesses['minOFF'] ])
        potential_fitnesses['mismatch']  = mismatch_score  

        t_sppms    = coarse_grain_pheno(bppms, base_to_seg, er_seg_ids, truth_vector, bool_TT)
        results.append( ( potential_fitnesses, t_sppms, mfe_dbs) )


    #print "Finished pheno perf fits"
    return results
