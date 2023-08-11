import numpy as np
np.set_printoptions(precision=4)
import math
import copy
from utility import*
import sys
import os
import pickle

class EA():
    def __init__(self):
        self.population = list()
        self.parents = list()
        self.offspring = list()
        self.id_counter = 0 #The 1st individual has id = 1 so that it is consistent with flk SQL primary key that starts indexing at 1          
               
    def initialize_population(self, initial_pop):
        if initial_pop == None:
            for i in range(self.num_inds):
                new_individual = self.create_individual()
                self.id_counter += 1 
                new_individual.id = self.id_counter 
                self.offspring.append(new_individual)
        else:
            for ind in initial_pop:
                self.offspring.append(self.copy_ind(ind))
            self.id_counter = len(initial_pop)

    def generate_semantics_offspring(self):
        pass

    def evaluate_fitness_offspring(self):
        pool = self.offspring
        #Calculate all potential fitnesses
        for current_individual in pool:
            current_individual.fitnesses = list()
            self.evaluate_fitness(current_individual)

        for ind in pool:
            ind.fitnesses = list()
            for fit_name in self.selected_fitnesses:
                ind.fitnesses.append(ind.potential_fitnesses[fit_name])

    def select_parents(self):
        if (self.parent_selection_method == 0): #Uniform parent selection
            self.parents = [self.copy_ind(ind) for ind in self.population]
        if (self.parent_selection_method == 2): #Tournament selection
            self.parents = [self.copy_ind(ind) for ind in self.population]
            self.parents = self.tournament_selection(self.parents)

    def tournament_selection(self,pool):
        winners = list()
        for i in range(self.num_inds): 
            tournament_winner = self.perform_tournament_round(pool)
            winners.append(self.copy_ind(tournament_winner))
        return winners

    def perform_tournament_round(self,pool):
        participants = list()
        for i in range(self.num_tournament_participants):
            current_participant_id = np.random.randint(0,self.num_inds)
            #This is with replacement, which reduces selection pressure e.g. if I have inds, the strongest one won't quite saturate immediately
            current_participant = pool[current_participant_id] 
            participants.append(current_participant)

        
        dataset = [p.fitnesses for p in participants]
        fittest_front = non_dominated_sort(dataset)[0]
        winner = participants[fittest_front[0]] #arbitrarily take the first element of best front
        return winner 

    def mutate_parent(self,parent):
        pass
        
    def crossover_parents(self,parent_1,parent_2):
        pass

    def copy_ind(self, ind):
        pass

    def produce_offspring(self):
        self.offspring = [self.copy_ind(ind) for ind in self.parents]

        num_parents = len(self.parents)
        num_parents_pairs = num_parents // 2
        #We randomize the order of the population so that we perform crossover with a random individual
        ids = np.random.choice(num_parents,  num_parents, replace = False) 
        #print ids
        for i in range(num_parents_pairs):
            sample = np.random.random_sample()
            parent_1 = self.offspring[ids[2*i]]
            parent_2 = self.offspring[ids[2*i+1]]
            self.id_counter += 1
            parent_1.id = self.id_counter
            self.id_counter += 1
            parent_2.id = self.id_counter
            
            if (sample < self.crossover_prob):
                self.crossover_parents(parent_1,parent_2)
            else:
                self.mutate_parent(parent_1)
                self.mutate_parent(parent_2)

  
    def select_survivors(self):
        population_plus_offspring = self.population  + self.offspring
        if self.survivor_selection_method == 0: ##GENITOR: replace worst
            for ind in population_plus_offspring:
                ind.fitnesses = list()
                for sF in self.selected_fitnesses:
                    ind.fitnesses.append(ind.potential_fitnesses[sF])
                
            dataset                              = [p.fitnesses for p in population_plus_offspring]
            population_plus_offspring_fronts_ids = flatten(non_dominated_sort(dataset))[0:self.num_inds]
            self.population                      = [population_plus_offspring[id] for id in population_plus_offspring_fronts_ids]

    def post_process(self):
        pass


    def save_checkpoint(self):
        checkpoint_dir = "checkpoints"
        #if (self.cur_gen + 1) % self.checkpoint_rate == 0 or self.cur_gen == 0:
        if (self.cur_gen + 1) % self.checkpoint_rate == 0:

            if not os.path.exists(checkpoint_dir):
                os.makedirs(checkpoint_dir)

            checkpoint_file = os.path.join(checkpoint_dir, f'checkpoint_gen_{self.cur_gen}.pkl')

            # Delete any existing checkpoint files
            for f in os.listdir(checkpoint_dir):
                if f.endswith('.pkl'):
                    os.remove(os.path.join(checkpoint_dir, f))


            with open(checkpoint_file, 'wb') as f:
                pickle.dump(self.population, f)

    def load_checkpoint(self, checkpoint_dir = "checkpoints"):      

        if not os.path.exists(checkpoint_dir):
            return None, None

        checkpoints = [f for f in os.listdir(checkpoint_dir) if f.startswith('checkpoint_gen_')]
        if not checkpoints:
            return None, None

        # Find the latest checkpoint file
        checkpoint_files = [f for f in os.listdir(checkpoint_dir) if f.endswith('.pkl')]
        if not checkpoint_files:
            raise FileNotFoundError('No checkpoint files found')

        latest_checkpoint = max(checkpoint_files, key=lambda f: int(f.split('_')[2].split('.')[0]))


        with open(os.path.join(checkpoint_dir, latest_checkpoint), 'rb') as f:
            population = pickle.load(f)
            #The population was saved at the end of this generation
            generation = int(latest_checkpoint.split('_')[-1].split('.')[0]) 

        return population, generation



    def print_fitnesses(self,pool):
        fronts_indices = non_dominated_sort([p.fitnesses for p in pool])
        fittest_front = [pool[f] for f in fronts_indices[0]]
        for obj_id in range(len(pool[0].fitnesses)):
            print("Objective " + str(obj_id))
            print([round(ind.fitnesses[obj_id],2) for ind in fittest_front])



    def check_existing_checkpoint_and_prompt(self, checkpoint_dir='checkpoints'):
        # Check if the checkpoint directory exists and contains checkpoint files
        if os.path.exists(checkpoint_dir) and any(f.startswith('checkpoint_gen_') and f.endswith('.pkl') for f in os.listdir(checkpoint_dir)):
            print('A checkpoint already exists. Do you want to:')
            print('  1. Continue from the saved checkpoint')
            print('  2. Start a new run')
            choice = input('Enter your choice (1 or 2): ').strip()

            if choice == '1':
                generation, population = self.load_checkpoint(checkpoint_dir)
                return generation, population
            elif choice == '2':
                # Delete all checkpoint files
                for f in os.listdir(checkpoint_dir):
                    if f.startswith('checkpoint_gen_') and f.endswith('.pkl'):
                        os.remove(os.path.join(checkpoint_dir, f))

                return None, None  # Start from generation 0
            else:
                print('Invalid choice. Exiting...')
                sys.exit(1)
        else:
            return None, None  # Start from generation 0

    def print_fitnesses(self,pool):
        ######
        fronts_indices = non_dominated_sort([p.fitnesses for p in pool])
        fittest_front = [pool[f] for f in fronts_indices[0]]
        for obj_id in range(len(pool[0].fitnesses)):
            print("Objective " + str(obj_id))
            print([round(ind.fitnesses[obj_id],2) for ind in fittest_front])



    def run(self, num_inds, num_gens, parent_selection_method, parent_selection_params, survivor_selection_method, crossover_prob, mutation_rate, mutation_probs, selected_fitnesses, bool_store):
        self.bool_store                = bool_store
        self.checkpoint_rate           = 3 #10 #1
        self.selected_fitnesses        = selected_fitnesses
        self.num_inds                  = num_inds
        self.num_gens                  = num_gens
        self.survivor_selection_method = survivor_selection_method #0 is genitor
        self.parent_selection_method   = parent_selection_method
        self.crossover_prob            = crossover_prob
        self.mutation_rate             = mutation_rate
        self.mutation_probs            = mutation_probs
        self.bool_store                = bool_store
        self.cur_gen                   = 0

        if self.parent_selection_method == 2:
            self.num_tournament_participants = int(parent_selection_params[0]) #in this case selection_params is a list of a single number

        #initial_pop, last_completed_gen = self.load_checkpoint()
        initial_pop, last_completed_gen = self.check_existing_checkpoint_and_prompt()

        if last_completed_gen == None:
            print("Starting generation " + str(self.cur_gen + 1) + " of " + str(self.num_gens))
            #print("Generation " + str(self.cur_gen))
            #print("Initializing population")
            #Treat the initial population like it is an offspring population. Initialize population creates an intial offspring but leaves population blank
            self.initialize_population(initial_pop) 
            #print("Generating semantics")
            self.generate_semantics_offspring()
            #print("Evaluating fitnesses") 
            self.evaluate_fitness_offspring()
            #print("Selecting survivors")
            self.select_survivors() #In this case, they are all survivors
            #self.print_fitnesses(self.population)
            if self.bool_store == 1:
                self.save_checkpoint()
            last_completed_gen = 0
        else:
            self.population = [self.copy_ind(ind) for ind in initial_pop]

        for current_generation in range(last_completed_gen+1, self.num_gens):
            self.cur_gen = current_generation
            print("Starting generation " + str(self.cur_gen + 1) + " of " + str(self.num_gens))
            #print("Generation " + str(current_generation))
            #print("Selecting parents")
            self.select_parents()
            #print("Producing offspring")
            self.produce_offspring()
            #print("Generating semantics")
            self.generate_semantics_offspring()
            #print("Evaluating fitnesses")
            self.evaluate_fitness_offspring()
            #print("Selecting survivors")
            self.select_survivors()
            if current_generation == self.num_gens - 1:
                self.post_process()

            #self.print_fitnesses(self.population)
            if self.bool_store == 1:   
                self.save_checkpoint()
        print("Finished generating selective ribozymes")