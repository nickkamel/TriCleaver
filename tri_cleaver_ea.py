from ea import *
from repr import SRS, SCC
from utility import multicore, search_list_objs, non_dominated_sort, space_out_items
from semantics import *
from math import ceil
import csv

class TriCleaver(EA):
    def __init__(self, viability_params, bool_fast_folding, folding_temp, downstream_seq, num_repeats_wt, num_repeats_mutant, repeat_unit, template_params):
        EA.__init__(self)       
        self.mutation_types = ['BBP', 'SBS']  #!!! changes this to: self.mutation_types = ['BCC', 'RzBS']     

        self.bool_fast_folding   = bool_fast_folding
        self.folding_temp        = folding_temp
        self.bool_TT             = 1 #Hardcoded

        self.min_cap, self.max_cap, self.bool_knee, self.knee_cap, self.knee_gen = viability_params

        self.spacing, self.spacing_length, self.s1_len, self.s3_len, self.sbs_len, max_num_strides, allowed_cut_site_types = template_params
        self.s2_len = 7 #!!! Hardcoded. TO DO: Specify this in template params Server side
        self.linker_len = 7 #!!! Hardcoded
        self.rzbs_len = self.s1_len + 1 + self.s3_len
        self.rzbs_extension_len = 10 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        self.misc_lens = self.s1_len, self.s2_len, self.s3_len, self.linker_len, self.rzbs_len, self.rzbs_extension_len


        self.wt_len = num_repeats_wt*3
        self.stride_size = ceil( (num_repeats_mutant - num_repeats_wt) / float(max_num_strides)    )*3
        self.repeats = repeat_unit*num_repeats_mutant

        #Build the list of allowed substrate binding sites
        #rs is the index of the cut site w.r.t. to the start of the sequence downstream of the repeats
        self.downstream_seq = downstream_seq
        self.allowed_rs = list()
        for i in range( len(self.downstream_seq) - self.sbs_len ):
            cut_site = self.downstream_seq[i + self.s3_len - 2:i + self.s3_len + 1]
            if cut_site in allowed_cut_site_types:
            #if self.downstream_seq[i + self.s3_len - 1] == 'U' and self.downstream_seq[i + self.s3_len] != 'G': #NUH
                self.allowed_rs.append(i)

        #Hardcode mismatch locations on sensor to give ribozyme some flexiblity when bound transcript 
        sensor_len = len(self.repeats)
        self.sensor_mismatches = list()
        for i in range(self.spacing_length):
            self.sensor_mismatches += space_out_items(i, self.spacing, sensor_len)

    def create_individual(self):
        ind                = SRS()
        ind.rs             = np.random.choice(self.allowed_rs)
        ind.update_representation(self.downstream_seq, self.sensor_mismatches, self.repeats, self.misc_lens, 1)
   
        return ind

    def copy_ind(self, ind):
        new_ind = copy.deepcopy(ind)
        return new_ind

    def pre_insert(self, inds):
        new_inds = list()
        for ind in inds:
            new_ind = ind
            new_ind.folding_constraints     = [task[1] for task in ind.tasks[0]]
            new_inds.append(new_ind)
        return new_inds

    def generate_semantics_offspring(self):
        population_full = [search_list_objs(self.population, 'id', ind.id)[0] for ind in self.parents] #Because I removed the deep copies, the parents don't have any phenotype

        if self.bool_fast_folding == 0 or (self.bool_fast_folding == 1 and self.cur_gen == self.num_gens - 1):
            bool_300 = 1
            pool    = population_full + self.offspring #We need to make sure ALL individuals in the final population go through bool 300 folding
        else:
            bool_300 = 0
            pool    = self.offspring

        batch = list()
        for ind in pool:
            #Prepare tasks
            ind.generate_segments()
            ind.generate_strands()
            ind.update_base_maps()
            ind.build_color_string()
            ind.generate_folding_tasks(self.wt_len, self.stride_size)
            batch.append(ind.tasks) 

        #Since all individuals have the same template, base_to_seg, name_to_bases, and er_seg_ids are the same for all individuals so we just take the ones from the last individual. Also note that unlike TS we do not specify a truth vector. Rather, the truth vector is determined by the number wt_len, sensor_len and num_strides. For now, the truth vector is determined for individual when the folding tasks are generated. Therefore TC has ind.truth_vector instead of self.truth_vector
        common_data = (ind.base_to_seg, ind.name_to_bases, ind.er_seg_ids, ind.truth_vector, bool_300, self.folding_temp, self.bool_TT, ind.sccs[4].mismatches)

        num_processes              = 12 #16
        self.num_processes_novelty = num_processes
        results                    = multicore(batch, pheno_perf_fits, num_processes, common_data)
        #results                   = pheno_perf_fits([batch, common_data]) #Single core


        #Store results back in individuals and prepare distance calculation
        phenos = list()
        for i, ind in enumerate(pool):
            result = results[i]
            ind.potential_fitnesses = result[0]
            ind.coarse_pheno        = result[1]
            ind.mfe_dbs             = result[2]

        #Novelty
        print("Starting novelty")        
        novelty_pool    = population_full + self.offspring
        num_inds        = len(novelty_pool)
        distance_matrix = self.calculate_distance_matrix(novelty_pool)  

        neighborhood_size = 30
        for i in range(num_inds):
            distances = list()
            for j in range(num_inds):
                distances.append(distance_matrix[i, j])                      

            neighbor_idxs                                  = arg_sort(distances)[0:neighborhood_size]
            novelty_pool[i].potential_fitnesses['novelty'] = sum([distances[n_idx] for n_idx in neighbor_idxs])


        #Viability nullification
        if self.bool_knee == 1:
            if self.cur_gen < self.knee_gen:
                self.cur_viability_threshold = self.min_cap + self.cur_gen * (self.knee_cap - self.min_cap) / float(self.knee_gen)
            else:
                self.cur_viability_threshold = self.knee_cap + (self.cur_gen - self.knee_gen) * (self.max_cap - self.knee_cap) / float(self.num_gens -1 - self.knee_gen)
        else:
            self.cur_viability_threshold = self.min_cap + self.cur_gen * (self.max_cap - self.min_cap) / float(self.num_gens -1)
        #print(self.cur_viability_threshold)     
        
        pool = self.population + self.offspring
        pfNames = list(ind.potential_fitnesses.keys())
        for ind in pool:
            for pfName in pfNames:
                if 'Cut' not in pfName: #Without this condition, we end up with nameCutCutCut ... Cut
                    if ind.potential_fitnesses['switchMin'] < self.cur_viability_threshold:
                        ind.potential_fitnesses[pfName + "Cut"] = -1000
                    else:
                        ind.potential_fitnesses[pfName + "Cut"] = ind.potential_fitnesses[pfName]

        
    def calculate_distance_matrix(self, pool):
        batch   = list()
        num_inds = len(pool)
        for i in range(num_inds): 
            for j in range(num_inds):
                batch.append((pool[i].coarse_pheno, pool[j].coarse_pheno))
        results = multicore(batch, dot_distance_batch, self.num_processes_novelty, [])
        #results = dot_distance_batch(batch) #Single core

        task_id               = 0
        distance_matrix = np.zeros(( num_inds, num_inds )) #Used later for MDS
        for i in range(num_inds):
            for j in range(num_inds):
                distance_matrix[i, j] = results[task_id]
                task_id += 1   

        return distance_matrix

    def post_process(self):
        #Acceptable
        acceptable_pool = list()
        for ind in self.population:
            #if ind.potential_fitnesses['switchMin'] >= -1: #!!!! Debugging !!!!!!
            if ind.potential_fitnesses['switchMin'] >= 0.85:
                ind.potential_fitnesses['acceptable'] = 1
                acceptable_pool.append(ind)
            else:
                ind.potential_fitnesses['acceptable'] = 0

        #Marginal diversity gain(MDG)
        #Handle unacceptable individuals
        for ind in self.population:
            if ind.potential_fitnesses['acceptable'] == 0:               
                ind.potential_fitnesses['mdg']        = -len(acceptable_pool) #This is effectively -inf

        #Handle acceptable individuals
        pool              = acceptable_pool
        num_inds          = len(pool)
        if num_inds == 0: #no acceptable individuals
            pass 
        elif num_inds == 1:
            pool[0].potential_fitnesses['mdg'] = 0
        else:
            distance_matrix   = self.calculate_distance_matrix(pool)      
            mdg_idxs          = list()
            furthest_pair_idx = np.argmax(distance_matrix)                          #This is the (single) index of the pair in a 2d matrix
            mdg_idxs += [furthest_pair_idx % num_inds, furthest_pair_idx // num_inds] #From this pair index, we get the index of each element of the pair
            while len(mdg_idxs) < num_inds:
                candidate_distances = list()
                unplaced_idxs       = [idx for idx in range(num_inds) if idx not in mdg_idxs]
                for unplaced_idx in unplaced_idxs:
                    marginal_distances = list()
                    for mdg_idx in mdg_idxs:
                        marginal_distances.append(distance_matrix[unplaced_idx, mdg_idx])
                    candidate_distances.append( min(marginal_distances)  )
                temp_idx = arg_sort(candidate_distances)[-1]
                mdg_idxs.append( unplaced_idxs[temp_idx] )

            pool[ mdg_idxs[0] ].potential_fitnesses['mdg'] = 0
            pool[ mdg_idxs[1] ].potential_fitnesses['mdg'] = 0
            for mdg_rank in range(2, num_inds):
                pool[ mdg_idxs[mdg_rank] ].potential_fitnesses['mdg'] = -mdg_rank

        #Export to CSV
        mdg_vals        = [ind.potential_fitnesses['mdg'] for ind in pool]
        ind_idxs_sorted = arg_sort(mdg_vals)[::-1]
        inds_ranked     = [pool[idx] for idx in ind_idxs_sorted]

        designs_folder = 'designs'
        if not os.path.exists(designs_folder):
            os.makedirs(designs_folder)

        header = ['Design ID', 'SRz sequence']
        rows = list()
        for ind_idx, ind in enumerate(inds_ranked):
            rows.append([ind_idx] + [ind.strands[0]])

        with open('designs/Srz_designs_.csv','w') as out:
            csv_out=csv.writer(out)
            csv_out.writerow(header)
            csv_out.writerows(rows)
           
    def evaluate_fitness(self, ind):
        pass

    def mutate_parent(self, ind):
        for m in range(self.mutation_rate):
            mutation_type = np.random.choice(self.mutation_types, 1, p=self.mutation_probs)[0]
            if mutation_type == 'SBS':
                options   = [i for i in self.allowed_rs if i != ind.rs]
                ind.rs    = np.random.choice(options)
                ind.update_representation(self.downstream_seq, self.sensor_mismatches, self.repeats, self.misc_lens, 1)
            else:
                #Randomly select a (global) index of a bcc. Use the bcc_to_scc map to get the scc containing that bcc and the bcc's relative position within that scc
                num_bccs_total       = len (ind.bcc_to_scc.keys() )
                bcc_id_all           = np.random.choice( range( num_bccs_total),  1, p=ind.mutation_prob)[0]
                (scc_id, bcc_id_scc) = ind.bcc_to_scc[bcc_id_all] 
                scc                  = ind.sccs[scc_id]
                bcc                  = scc.bccs[bcc_id_scc]
                options              = [i for i in scc.valid_nucs[bcc_id_scc] if i != bcc] #Every bcc but the current one      
                scc.bccs[bcc_id_scc] = np.random.choice(options)     

                
