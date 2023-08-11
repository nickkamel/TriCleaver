from multiprocessing import Pool
import numpy as np

def multicore(tasks, func, num_processes, common_data):
    num_tasks_per_process_floor = len(tasks) // num_processes #This is rounded to nearest integer
    num_tasks_total_floor       = num_tasks_per_process_floor*num_processes
    tasks_whole                 = tasks[0:num_tasks_total_floor]
    tasks_remainder             = tasks[num_tasks_total_floor:len(tasks)]

    results = list()
    if len(tasks_whole) !=0:
        n = num_tasks_per_process_floor
        #processes_tasks = [tasks_whole[i:i + n] for i in range(0, num_tasks_total_floor, n)]
        processes_tasks = [(tasks_whole[i:i + n], common_data) for i in range(0, num_tasks_total_floor, n)]
        pool = Pool(num_processes) 
        processes_results = pool.map(func, processes_tasks)
        for p in processes_results:
            results += p
        pool.close()
        pool.join()

    if len(tasks_remainder) != 0:
        results += func([tasks_remainder, common_data])

    return results

def prefix_sum(input):
    accum = 0
    out = list()
    for i in input:
        accum += i
        out.append(accum)
    return out

def avg(input_list):
    return sum(input_list)/float(len(input_list))

#Sorts in ascending order 
def arg_sort(seq):
    return [i for (v, i) in sorted((v, i) for (i, v) in enumerate(seq))]

def flatten(l):
    flat_list = list()
    for sublist in l:
        for item in sublist:
            flat_list.append(item)
    return flat_list

def search_list_objs(target_list, attribute, target_value):
    results = list()
    for stored_object in target_list:
        if getattr(stored_object, attribute) == target_value:
            results.append(stored_object)
    return results 

def inv_list(target_list, value):
    results = [i for i, x in enumerate(target_list) if x == value]
    if results != []:
        return results[0]
    else:
        return None

#Takes a list of lists, where the lower level list represents the values for the multiple objectives
def non_dominated_sort(dataset):
    #print dataset
    num_attributes = len(dataset[0])
    list_size = len(dataset)
    num_ranked_individuals = 0
    dominated_counts = list()
    domination_sets = list()
    fronts = list()
    front0 = list()
    for i in range(list_size):
        domination_set = list()
        dominated_count = 0
        for j in range(list_size):
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
            #print dataset[i], dataset[j], num_better_attributes_i, num_better_attributes_j, num_equal_attributes
            
            if num_better_attributes_j == 0 and num_equal_attributes != num_attributes: #if j is never better than i and i is better than j at at least 1 fitness
                 domination_set.append(j)
            elif num_better_attributes_i == 0 and num_equal_attributes != num_attributes: #j dominates i
                 dominated_count += 1
                        
        if dominated_count == 0:
            front0.append(i)
            num_ranked_individuals += 1
        domination_sets.append(domination_set)
        dominated_counts.append(dominated_count)

    fronts.append(front0)
    front_counter = 0
    while(num_ranked_individuals != list_size):
        new_front = list()
        for i in fronts[front_counter]:
            for j in domination_sets[i]:
                dominated_counts[j] -= 1
                if dominated_counts[j] == 0:
                    new_front.append(j)
                    num_ranked_individuals +=1                        
        fronts.append(new_front)
        front_counter +=1

    return fronts

def dot_distance(v1, v2):
    return np.sum(np.abs(np.subtract(v1, v2)) )

def dot_distance_batch(batch):
    tasks, common_data = batch
    results = list()
    for task in tasks:
        results.append(dot_distance(*task))
    return results

def unique(target_list):
    unique_list = list(set(target_list)) #Doesn't maintain order of list
    return unique_list

def bool_func_int_to_list(bool_func_int, num_bits):
    format_string = '{0:0' + str(num_bits) + 'b}'
    #temp = list('{0:08b}'.format(bool_func_int))
    temp = list(format_string.format(bool_func_int))
    boolFuncList = [int(c) for c in temp]
    return boolFuncList

def bin_to_dec(binary):
    return sum(val*(2**idx) for idx, val in enumerate(reversed(binary)))
    
def Remove_From_List(input_list, items):
    for item in items:
        if item in input_list:
            input_list.remove(item)
    return input_list

def space_out_items(offset, spacing, list_len):
    item_idxs = list()
    item_idx  = offset + spacing
    while item_idx < list_len:
        item_idxs.append(item_idx)
        item_idx = item_idx + 1 + spacing
    return item_idxs

def convert_to_float(s):
    try:
        float_value = float(s)
        return float_value
    except ValueError:
        return None