from utility import avg

#This functions assigns a score from 0 (no disruption) to 1 (full disruption) to each motif
#In this case for TC, the motifs are C1, C2, C3, stem1, stem2, stem3. In TS C1, C2, C3 are treated as a single core motif
def evaluate_rz_activity(bppm, name_to_bases):  
    #core
    core_scores = list()
    core_score_abs = 0
    for core_seg_name in ['C1', 'C2', 'C3']:
        core_seg_score_abs = 0
        core_seg_base_ids  = name_to_bases[core_seg_name]
        for base_id in core_seg_base_ids:
            unpaired_prob = bppm[(base_id, None)]
            paired_prob   = 1 - unpaired_prob
            core_seg_score_abs += paired_prob
        core_seg_score = core_seg_score_abs / float(len(core_seg_base_ids))
        core_scores.append(core_seg_score)
        core_score_abs += core_seg_score_abs
 
    num_core_base_ids  = len(name_to_bases['C1'] + name_to_bases['C2'] + name_to_bases['C3'])
    core_score         = core_score_abs / num_core_base_ids #This will be slightly different than the average of core_scores

    #stems
    stem_scores = list()
    for i in range(1, 4):
        stem_name     = 's' + str(i)
        stem_base_ids = list(zip(name_to_bases[stem_name + 'A'], name_to_bases[stem_name + 'B'][::-1]))
        stem_score_  = 0
        stem_score_A = 0
        stem_score_B = 0
        w_up = 1 #unpaired weight
        w_mp = 1 #paired weight
        for baseId, partner_base_id in stem_base_ids: 
            paired_prob             = bppm[(baseId, partner_base_id)]
            mispaired_unpaired_prob = 1 - paired_prob
            unpaired_prob1          = bppm[(baseId, None)]
            mispaired_prob1         = mispaired_unpaired_prob - unpaired_prob1
            unpaired_prob2          = bppm[(partner_base_id, None)]
            mispaired_prob2         = mispaired_unpaired_prob - unpaired_prob2
            stem_score_A            = stem_score_A + w_up*unpaired_prob1 + w_mp*mispaired_prob1
            stem_score_B            = stem_score_B + w_up*unpaired_prob2 + w_mp*mispaired_prob2

        stem_score_A = stem_score_A / float( len(stem_base_ids) )
        stem_score_B = stem_score_B / float( len(stem_base_ids) )
        stem_score  = max(stem_score_A, stem_score_B)
        #stem_score  = avg([stem_score_A, stem_score_B]) #TS
        stem_scores.append(stem_score)

    #rz_elem_scores = [core_score] + stem_scores #TS
    rz_elem_scores = core_scores + stem_scores

    on_score  = -avg(rz_elem_scores) #This reaches max value when the core and stem motifs experience no disruption
    off_score = max(stem_scores)     #This reaches max value when the at least one stem is completely disrupted

    return on_score, off_score