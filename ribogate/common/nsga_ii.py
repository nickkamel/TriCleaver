def non_dominated_sort(num_items, num_attributes, dataset):
    '''
    Implements the NSGA-ii multiobjective sorting algorithm from Deb, Kalyanmoy, et al. "A fast and elitist multiobjective genetic algorithm: NSGA-II."
    Sorts the items in a series of non-dominated fronts
    Each item in the same front is dominated by the same number of items
    An item A dominates an item B if 1) A is better than B at least one objective 2) A is no worse than B at any objective
    Examples: (5,5) dominates (5,0). Neithher (5,0) and (0,5) dominate each other.
    
    Parameters:
        num_items (int): The number of items to sort
        num_attributes (int): The number of attributes (e.g. objectives) that each item contains. If it only contains one, NSGA-ii is equivalent to simply sorting a list of numbers
        dataset (num_items lists each containing num_attribtes ints/floats)
    
    Returns 
        fronts (list of lists of int). Each element in the higher level list corresponds to a non-dominated front, which in turn is a list contained the indices of the items that are part of that front
    '''

    num_ranked_items = 0
    dominated_counts = list()
    domination_sets = list()
    fronts = list()
    top_front = list()
    for i in range(num_items):
        domination_set = list()
        dominated_count = 0
        for j in range(num_items):
            num_better_attributes_i = 0
            num_better_attributes_j = 0
            num_equal_attributes = 0
            for attributeId in range(num_attributes):
                if dataset[i][attributeId] > dataset[j][attributeId]:
                    num_better_attributes_i += 1
                elif dataset[i][attributeId] < dataset[j][attributeId]:
                    num_better_attributes_j += 1
                else:
                    num_equal_attributes += 1
            
            if num_better_attributes_j == 0 and num_equal_attributes != num_attributes: # if j is never better than i and i is better than j at at least 1 fitness
                 domination_set.append(j)
            elif num_better_attributes_i == 0 and num_equal_attributes != num_attributes: # j dominates i
                 dominated_count += 1
                        
        if dominated_count == 0:
            top_front.append(i)
            num_ranked_items += 1
        domination_sets.append(domination_set)
        dominated_counts.append(dominated_count)

    fronts.append(top_front)
    front_counter = 0
    while(num_ranked_items != num_items):
        new_front = list()
        for i in fronts[front_counter]:
            for j in domination_sets[i]:
                dominated_counts[j] -= 1
                if dominated_counts[j] == 0:
                    new_front.append(j)
                    num_ranked_items +=1                        
        fronts.append(new_front)
        front_counter +=1

    return fronts
