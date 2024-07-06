from multiprocessing import Pool

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

def arg_sort(input_list):
    '''
    Arg sorts in ascending order 
    '''
    
    return [i for (v, i) in sorted((v, i) for (i, v) in enumerate(input_list))]

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

def unique(target_list):
    unique_list = list(set(target_list)) #Doesn't maintain order of list
    return unique_list

def bool_func_int_to_list(bool_func_int, num_bits):
    format_string = '{0:0' + str(num_bits) + 'b}'
    #temp = list('{0:08b}'.format(bool_func_int))
    temp = list(format_string.format(bool_func_int))
    bool_func_list = [int(c) for c in temp]
    return bool_func_list

def bin_to_dec(binary):
    return sum(val*(2**idx) for idx, val in enumerate(reversed(binary)))
    
def remove_items_from_list(input_list, items):
    '''
    inputs:
        input_list (list): A list from which we want to remove multiple items
        items (list): A list of items that we want to remove from the input_List
    return:
        input_list (list): The original input list with the items removed
    '''

    for item in items:
        if item in input_list:
            input_list.remove(item)
    return input_list

def convert_to_float(s):
    try:
        float_value = float(s)
        return float_value
    except ValueError:
        return None