import numpy as np
import copy

import config
from ..common.ea import EA
from .ribozyme_templates import RibozymeTemplate
from ..common.representation import SwitchableRibozyme, initialize_bccs
from .representation import create_representation, update_representation
from .folding_preparation import generate_truth_vector, calculate_stride_size
from ..common.mutation import mutate_bcc
from ..common.novelty import calculate_novelty_scores
from .blast import run_blast, calculate_specificity
from .fitness_evaluation import predict_structures_evaluate_performance, fold_unhybridized_srz, calculate_distance_matrix_tc
from ..common.utility import multicore
from ..common.reporting import determine_reported_individuals_and_marginal_diversity_rank

class TriCleaver(EA):
    def __init__(self, ribozyme_type, template_params, folding_temperature, downstream_seq, num_repeats_wt, num_repeats_mutant, repeat_unit, allowed_rzbs_positions): 
        '''
        Initialize for the TriCleaver class which inherits from the EA class and overrides several of its functions.
        Configures TriCleaver for the user-specified trinucleotide(*) repeat expansion disorder (TRED)
        (*) It can now also handle diseases with repeat lengths other than 3, although more testing should be done for these cases
        Note that extraction of valid ribozyme binding sites (RzBSs) has already been performed when validating user input
        Therefore, this function's arguments are not the raw WT and mutant sequences, but data that has been extracted from them 

        Parameters:
            ribozyme_type (string): The name of the ribozyme that we are using. Currently two options: 'minimal_hammerhead' and 'extended_hammerhead'
            template_params (dict): Contains parameters relating to the ribozyme's structure  
            folding_temperature (float)
            downstream_seq (string): The substrate sequence downstream of the repeats segment. Identical between the wild-type and mutant
            num_repeats_wt (int): Number of repeats on the wild-type strand
            num_repeats_mutant (int): Number of repeats on the mutant strand. See after parameters list
            repeat_unit (string): The seqeuence that is repeated (e.g. CAG)
            allowed_rzbs_positions (list of ints): Positions w.r.t. start of downstream_seq of potential RzBS with cut site motifs okayed by the user

        The mutant repeats stretch considered by TriCleaver can differ slightly from the actual repeats segment on the substrte for three possible reasons:
        1) We cap the number of repeats if they are very long
        2) If the repeats consist of mutiple codons mapping to the same amino acid, we only consider the most common codon
        3) TriCleaver is not bound by the open reading frame. 
        Consider the case of 4 (GCG) repeats, flanked by a G and GC, respectively: G-(GCG)(GCG)(GCG)(GCG)-GC
        When TriCleaver scans for repeats, if will actually detect 5 GGC repeats (GGC)(GGC)(GGC)(GGC)(GGC)  

        '''
        
        EA.__init__(self)       
        self.mutation_types = ['BCC', 'RzBS']      

        self.folding_temperature = folding_temperature
        self.are_segments_low_level = template_params['ARE_SEGMENTS_LOW_LEVEL']
        self.rzbs_diversity_weight = template_params['RZBS_DIVERSITY_WEIGHT']
        self.repeat_unit_len = len(repeat_unit)
        self.num_repeats_wt = num_repeats_wt
        self.num_repeats_mutant = num_repeats_mutant
        self.wt_len = self.num_repeats_wt*self.repeat_unit_len
        self.mutant_len = self.num_repeats_mutant*self.repeat_unit_len
        self.sensor_len = self.mutant_len
        self.repeats = repeat_unit*self.num_repeats_mutant 
        self.downstream_seq = downstream_seq
        self.allowed_rzbs_positions = allowed_rzbs_positions # These positions are defined w.r.t. to the start of the sequence downstream of the repeats

        # Load in the relevant template
        self.ribozyme_type = ribozyme_type
        self.template_params = template_params
        self.template = RibozymeTemplate(ribozyme_type, template_params, self.downstream_seq, self.repeats)

        self.stride_size = calculate_stride_size(self.num_repeats_wt, self.num_repeats_mutant, self.template.max_num_strides, 3)
        self.truth_vector = generate_truth_vector(self.sensor_len, self.wt_len, self.stride_size)
        self.sensor_mismatches = self.template.scc_mismatches[('SENSOR','REPEATS')]
        
        # This data is constant between individuals and throughout the EA and it is required for fitness evaluation
        self.constant_data = self.template.extension_region_seg_ids, self.truth_vector, self.folding_temperature, self.are_segments_low_level, self.sensor_mismatches, config.folding_package_paths, self.template.ribozyme_type, self.template, self.wt_len, self.stride_size
                 
    def create_individual(self):     
        # Instantiates a candidate switchable ribozyme individual
        ind = SwitchableRibozyme()
        
        # Chooses a random ribozyme binding site
        ind.rzbs_position = np.random.choice(self.allowed_rzbs_positions)
        
        # Creates a set of segment connected components (SCCs) based on that site and the switchable ribozyme template
        # Also creates a bcc_to_scc map and a list of mutation_probabilities
        ind.sccs, ind.bcc_to_scc, ind.bcc_mutation_probabilities = create_representation(self.template, ind.rzbs_position)
        
        # Samples random valid base connected components (bccs) for each SCC
        for scc in ind.sccs:
            scc.bccs = initialize_bccs(scc.length, scc.valid_nucs)
        return ind

    def copy_individual(self, ind):
        new_ind = copy.deepcopy(ind)
        return new_ind

    def pre_insert(self, inds):
        new_inds = list()
        for ind in inds:
            new_ind = ind
            new_inds.append(new_ind)
        return new_inds

    def evaluate_fitness(self):
        '''      
        Fitness evaluation is performed in two steps:
        1) The structures of each offspring are predicted in parallel. The perforamnce fitness measures (e,g. ON, OFF, mismatch), are also calculated in parallel.
        2) The results are returned and then novelty is calculated between the offspring and parents
        '''       
        # These two functions update the individuals in place 
        predict_structures_evaluate_performance(self.offspring, config.folding_packages_tc['LIVE'], self.constant_data)        
        calculate_novelty_scores(self.population + self.offspring, self.novelty_neighborhood_size, calculate_distance_matrix_tc, {'rzbs_diversity_weight':self.rzbs_diversity_weight})
            
    def mutate_individual(self, ind):    

        # Applies a mutation mutation_rate times
        for m in range(self.mutation_rate):
            
            # For each mutation, randomly selects a mutation type from the self.mutation_probabilities distribution     
            mutation_type = np.random.choice(self.mutation_types, 1, p=self.mutation_probabilities)[0]
            
            # Change to another valid ribozyme binding site
            # This requires us to update the representation 
            if mutation_type == 'RzBS':
                options = [i for i in self.allowed_rzbs_positions if i != ind.rzbs_position]
                ind.rzbs_position = np.random.choice(options)
                ind.sccs, ind.bcc_to_scc, ind.bcc_mutation_probabilities = update_representation(self.template, ind.rzbs_position, ind.sccs, False) 
                
            # Mutate a random bcc
            elif mutation_type == 'BCC':
                mutate_bcc(ind.bcc_to_scc, ind.bcc_mutation_probabilities, ind.sccs)
                
            else:
                raise NotImplementedError(f"There is no mutation function for mutation type: ({mutation_type})")

                
    def post_process(self):
        '''
        This function performs a variety of functions that are either too costly to perform each generation during the EA or that server no purpose during the EA loop itself
        Currently performs the following computation: 
            * Folds the final population with a different set of folding packages for more accurate prediction
            * Calculates a specificity score by blasting each individual against the human genome.
            * Folds each sRz on their own (i.e. without the substrate).
            * Determines which individuals are acceptable based on certain criteria.
            * Calculates a marginal diversity rank for each acceptable individual.
        '''
        
        if config.folding_packages_tc['POST'] is not None:
            predict_structures_evaluate_performance(self.population, config.folding_packages_tc['POST'], self.constant_data)
            
            '''
            In previous implementations, I would fold the offspring of the last generation and the parents of the previous generation and then performed a final round of survivor selection on this merged population to get the final population
            However, this was confusing. The current approach is simpler, but I can't calculate novelty since that would require comparing to the parents as well. 
            Therefore I just set the novelty score to 0 since it carries no informative value on its own anyway (it only makes sense w.r.t. other individuals)    
            '''
            for ind in self.population:
                ind.potential_fitnesses['novelty'] = 0

        # Calculate specificity score by blasting the RzBS against the human genome
        # Currently only uses a single process     
        if config.should_perform_blast:
            for ind in self.population:
                rzbs = self.downstream_seq[ind.rzbs_position:ind.rzbs_position + self.template.rzbs_len]
                blast_output = run_blast(rzbs, config.blast_db_path)
                ind.potential_fitnesses['specificity'] = calculate_specificity(blast_output)
                #print("specificity", ind.potential_fitnesses['specificity'], rzbs)

        # Kharma wanted me to fold the SRzs on their own. This is not used as part of the EA loop. It is only used to visualize the sRz on its own on the web app. 
        common_data = (self.folding_temperature, config.folding_packages_tc['POST'], config.folding_package_paths)
        batch = [[ind.strands[0], '.'*len(ind.strands[0])] for ind in self.population]
        results = multicore(batch, fold_unhybridized_srz, config.num_processes, common_data)
        for idx, result in enumerate(results):
            self.population[idx].mfe_db_unhybridized_srz = result[0]
            
               
        # Fold the sRZs with a longer region flanking the RzBS
        if self.template.rzbs_extension_lens['POST'] > self.template.rzbs_extension_lens['LIVE']:        
            #individuals = copy.deepcopy(self.population) # Use this if we don't want to actually change the members of the population, just get the fitness results with the longer region
            individuals = self.population
            for ind in individuals:
                ind.sccs, ind.bcc_to_scc, ind.bcc_mutation_probabilities = update_representation(self.template, ind.rzbs_position, ind.sccs, True)
                
            predict_structures_evaluate_performance(individuals, config.folding_packages_tc['POST'], self.constant_data)
            for ind in self.population:
                ind.potential_fitnesses['novelty'] = 0
                
        # Determine which individuals should be reported to the user. Rank these by marginal diversity.
        self.population = determine_reported_individuals_and_marginal_diversity_rank(self.population, config.reportability_criteria_tc, calculate_distance_matrix_tc, {'rzbs_diversity_weight':self.rzbs_diversity_weight})