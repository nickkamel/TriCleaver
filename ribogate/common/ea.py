import numpy as np
np.set_printoptions(precision=4)

from .utility import flatten
from .nsga_ii import non_dominated_sort


class EA():
    def __init__(self):
        self.population = list()
        self.parents = list()
        self.offspring = list()
        
        # The 1st individual has id = 1 so that it is consistent with flk SQL primary key that starts indexing at 1     
        self.id_counter = 0   
        
    def set_ea_params(self, ea_params, store_individuals_func, get_job_status_func, job_id):
        '''       
        Parameters:
            ea_params (dict): The general parameters used by the EA (agnostic to the specific use case of the EA)
                NUM_INDS (int): # of individuals in the population
                NUM_GENERATIONS (int)
                PARENT_SELECTION_METHOD (string)
                SURVIVOR_SELECTION_METHOD (string)
                CROSSOVER_PROBABILITY (float): Probabiltiy of choosing crossover vs mutation. Currently, always set to 0 (no crossover is performed)
                MUTATION_RATE (float): The extent to which an individual is mutated w.r.t. its parent
                MUTATION_PROBABILITIES (list of floats): The probability of choosing a specific mutation operator
                NOVELTY_NEIGHBORHOOD_SIZE (int): The # of nearest neighbors (in phenotype space) considered when calculating an individual's novelty score
                IS_VIABILITY_NULLIFICATION_ENABLED (bool): Whether or not to perform viability nullification
                VIABILITY_PARAMS (dict): The parameters used by viability nullification
                SELECTED_FITNESSES (list of strings): The names of the fitness terms that are actually used for parent and survivor selection
            store_individual_function (function): A function that saves the population to disk. 
                When running from web app, provide a function that inserts population into database
                When running from standalone CLI, provide a function that saves the population as pickled objects
                By passing the function this way, we can completely remove any database and web app imports here which faciliates installing the CLI version
            get_job_status_func (function): A function that retrieves the current status of the job
                Only relevant for the web app case. The web app runs in a seperate thread and if users cancel the jobs from that thread, ..
                    ... we need a way for the EA to retrieve this updated status
                For consistency, we just provide a dummy get_job_status_func for the CLI version
            job_id (int):
                Only needed for the web app case (to retrieve the job's status from the database)
                Pass a None value for the CLI case                
        '''
        self.num_inds = ea_params['NUM_INDS']
        self.num_generations = ea_params['NUM_GENERATIONS']
        self.parent_selection_method = ea_params['PARENT_SELECTION_METHOD']
        self.survivor_selection_method = ea_params['SURVIVOR_SELECTION_METHOD']
        self.crossover_prob = ea_params['CROSSOVER_PROBABILITY']
        self.mutation_rate = ea_params['MUTATION_RATE']
        self.mutation_probabilities = ea_params['MUTATION_PROBABILITIES']
        self.novelty_neighborhood_size = ea_params['NOVELTY_NEIGHBORHOOD_SIZE']       
        self.viability_params = ea_params['VIABILITY_PARAMS']
        
        
        # See comments in def apply_viability_nullification for why we appen 'Via' to fitness naems
        if ea_params['IS_VIABILITY_NULLIFICATION_ENABLED']:
            self.selected_fitnesses = list()
            for name in ea_params['SELECTED_FITNESSES']:
                self.selected_fitnesses.append(name + 'Via')
        else:
            self.selected_fitnesses = ea_params['SELECTED_FITNESSES']        

        # These parameters will be different depending on whether the EA is run from a web server or a standalone CLI
        self.store_individuals_func = store_individuals_func 
        self.get_job_status_func = get_job_status_func
        self.job_id = job_id

        # TO DO: If database operations become a bottleneck, set the insertion_rate to a higher value
        self.insertion_rate = 1

               
    def initialize_population(self, initial_pop):
        print("Initializing population")

        if initial_pop == None:
            for i in range(self.num_inds):
                new_individual = self.create_individual()
                self.id_counter += 1 
                new_individual.id = self.id_counter 
                self.offspring.append(new_individual)
        else:
            for ind in initial_pop:
                self.offspring.append(self.copy_individual(ind))
            self.id_counter = len(initial_pop)

    def evaluate_fitness(self):
        '''
        Evalutes the fitness of a group of individuals (usually the population). Implemented in the derived class
        '''
        pass
    
    def apply_viability_nullification(self):
        '''
        Raises the viability threshold as described in "Evolutionary design and analysis of ribozyme-based logic gates" by N.Kamel and thesis
        If an individual's viability is below this threshold, set all its fitnesses to "negative infinity"
        '''
        
        min_via, max_via = self.viability_params['MIN_VIA'], self.viability_params['MAX_VIA']

        final_gen = self.num_generations -1
        if self.viability_params['HAS_BREAKPOINT']:
            breakpoint_value, breakpoint_gen = self.viability_params['BREAKPOINT_VALUE'], self.viability_params['BREAKPOINT_GEN']
            if self.current_generation < breakpoint_gen:
                self.cur_viability_threshold = min_via + self.current_generation * (breakpoint_value - min_via) / float(breakpoint_gen)
            else:
                self.cur_viability_threshold = breakpoint_value + (self.current_generation - breakpoint_gen) * (max_via - breakpoint_value) / float(final_gen - breakpoint_gen)
        else:
            self.cur_viability_threshold = min_via + self.current_generation * (max_via - min_via) / float(final_gen)
            
        
        '''
        For each potential fitness term, we create a variant with "Via" appended to its name
        If the individual's viability is above the threshold, the Via variant is equal to the original fitness term
        Otherwise, set the Via variant to "negative infinity"
        This allows us to nullify fitnesses while still keep track of their original value for reporting
        '''
        individuals = self.population + self.offspring
        fit_names = list(individuals[0].potential_fitnesses.keys())
        for ind in individuals:
            for name in fit_names:
                if 'Via' not in name: # Without this condition, we end up with nameViaViaVia ... Via
                    if ind.potential_fitnesses['viability'] < self.cur_viability_threshold:
                        ind.potential_fitnesses[name + "Via"] = float('-inf') # Nullify fitness
                    else:
                        ind.potential_fitnesses[name + "Via"] = ind.potential_fitnesses[name]
                        
        num_viable_inds = len( [ind for ind in individuals if ind.potential_fitnesses['viability'] >= self.cur_viability_threshold] )
        print("Current viability threshold:", self.cur_viability_threshold, "# viable individuals:", num_viable_inds)     


    def apply_fitness(self):
        '''
        Copies of the values of certain selected fitnesses from the ind.potential_fitnesses dict to the ind.fitnesses list that is actually used for selection
        During fitness evaluation, multiple potential fitnesses are calculated
        However, not all of these will actually be used for selection during the EA
        Only those whose names are specified in selected_fitnesses are used
        '''

        for ind in self.offspring:
            ind.fitnesses = list()
            for fit_name in self.selected_fitnesses:
                ind.fitnesses.append(ind.potential_fitnesses[fit_name])

    def select_parents(self):
        '''
        Select parents from the population. These parents will later produce offspring.
        Currently implements one type of parent selection method: 'UNIFORM' where each individual in the population is selected as a parent
        '''
        print("Selecting parents")

        if self.parent_selection_method == 'UNIFORM':
            self.parents = [self.copy_individual(ind) for ind in self.population]
        else:
            raise NotImplementedError(f"There is no parent selection method called: ({self.parent_selection_method})")

    def mutate_individual(self, ind):
        '''
        Mutates an individual. Implemented in the derived class
        '''
        pass
        
    def crossover_individuals(self, ind_1, ind_2):
        raise NotImplementedError(f"Crossover is not implemented")

    def copy_individual(self, ind):
        '''
        For performance reasons, we may only want certain aspects (e.g. its genotype but not its phenotype) of an individual to be copied from generation
        Implemented in the derived class.
        '''
        pass

    def produce_offspring(self):
        '''
        Produces a set of offspring that are mutated copies of their parents
        '''
        print("Producing offspring")

        # First copy the parents and then mutate them (or perform crossover)
        self.offspring = [self.copy_individual(ind) for ind in self.parents]

        num_parents = len(self.parents)
        num_parents_pairs = num_parents // 2
        # We randomize the order of the population so that we perform crossover with a random individual
        ids = np.random.choice(num_parents,  num_parents, replace = False) 
        for i in range(num_parents_pairs):
            sample = np.random.random_sample()
            parent_1 = self.offspring[ids[2*i]]
            parent_2 = self.offspring[ids[2*i+1]]
            self.id_counter += 1
            parent_1.id = self.id_counter
            self.id_counter += 1
            parent_2.id = self.id_counter
            
            if (sample < self.crossover_prob):
                self.crossover_individuals(parent_1, parent_2)
            else:
                self.mutate_individual(parent_1)
                self.mutate_individual(parent_2)

  
    def select_survivors(self):
        '''
        Select a set of survivors from the parents and offspring to become the new population of the next generation
        '''
        print("Selecting survivors")

        # Merge the parents and offspring (they compete with each other to be potential survivors)
        population_plus_offspring = self.population  + self.offspring
        
        # Implement the selected survivor selection option
        # Currently there is one option: 'GENITOR', which culls the worst performing individuals
        if self.survivor_selection_method == 'GENITOR':
            for ind in population_plus_offspring:
                ind.fitnesses = list()
                for fit_name in self.selected_fitnesses:
                    ind.fitnesses.append(ind.potential_fitnesses[fit_name])                    
                 
            # Perform NSGA-ii ranking
            num_items, num_attributes = len(population_plus_offspring), len(self.selected_fitnesses)
            dataset = [p.fitnesses for p in population_plus_offspring]
            fronts_idxs = non_dominated_sort(num_items, num_attributes, dataset)
            
            # The survivors are the fittest half of the merged parent and offspring population
            fronts_idxs = flatten(fronts_idxs)[0:self.num_inds] 
            self.population = [population_plus_offspring[idx] for idx in fronts_idxs]           
        else:
            raise NotImplementedError(f"There is no survivor selection method called: ({self.survivor_selection_method})")

    def post_process(self):
        '''
        Performs additional processing on the set of individuals generated by the EA.
        This processing occurs once the EA has ended. Implemented in the derived class.
        '''
        pass


    def pre_insert(inds):
        '''
        To reduce storage requirements, we may want only to store certain attributes of the individuals.
        Implemented in derived class
        '''
        pass
    
    def store_individuals(self, should_overwrite_old_individuals, insertion_rate):
        self.store_individuals_func(self, should_overwrite_old_individuals, insertion_rate)

    def print_fitnesses(self, pool):
        '''
        Prints out the fitness values for each individual in the top front
        '''

        fronts_idxs = non_dominated_sort( len(pool), len(self.selected_fitnesses), [p.fitnesses for p in pool])
        top_front = [pool[f] for f in fronts_idxs[0]]
        for i in range(len(pool[0].fitnesses)):
            print("Objective " + str(i))
            print([round(ind.fitnesses[i],2) for ind in top_front])

    # No need to push an app context here since one is already pushed to the scheduler and this runs in the same thread as the scheduler
    # In the past I tried running this in a different thread than the scheduler and pushing its own app context but this messed up the database objects and trying to debug it was a disaster.
    def run(self, initial_pop, last_completed_generation):
        '''
        Runs the EA.   
        
        Parameters:
            iniital_pop (list of individuals): If initial_pop is not None, the EA skips initialization and resumes from this population
            last_completed_generation (int): The generation associated with initial_pop. -1 if initial_pop is None.      
        '''
        self.status = 'Running'
        self.current_generation = 0

        # If the EA is starting from scratch, we perform the iniital generation (which is slightly different from all subsequent generations due to initialization)
        if last_completed_generation == -1:
            print("Generation " + str(self.current_generation))
            #Treat the initial population like it is an offspring population. Initialize population creates an intial offspring but leaves population blank
            self.initialize_population(initial_pop)             
            print("Evaluating fitness") 
            self.evaluate_fitness()
            self.apply_viability_nullification()
            self.apply_fitness()            
            self.select_survivors() #In this case, they are all survivors            
            self.print_fitnesses(self.population)            
            
            self.store_individuals(True, self.insertion_rate)
            last_completed_generation = 0
        # If the EA is resuming from a previous state in the DB, then simply copy the provided initial population to the current population
        else:
            self.population = [self.copy_individual(ind) for ind in initial_pop]

        for current_generation in range(last_completed_generation + 1, self.num_generations):
            status = self.get_job_status_func(self.job_id)
            if status == 'Cancelled':
                break
            
            print("Generation " + str(current_generation))
            self.current_generation = current_generation
            
            self.select_parents()            
            self.produce_offspring()            
            print("Evaluating fitness") 
            self.evaluate_fitness()
            self.apply_viability_nullification()
            self.apply_fitness()            
            self.select_survivors()            
            print("Post processing")
            if current_generation == self.num_generations - 1:
                self.post_process()
            self.print_fitnesses(self.population)            
            self.store_individuals(True, self.insertion_rate)
                
        print("Finished EA")
        if status != 'Cancelled': 
            self.status = 'Completed'

        
