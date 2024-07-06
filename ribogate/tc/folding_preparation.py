from math import ceil
from ..common.folding_constraints import generate_folding_constraint_string

def generate_truth_vector(sensor_len, wt_len, stride_size):
    '''
    Builds the target output state (cleavage or no cleavage) for each possible hydrization state
    See supplementary material section 2.4 of "An Evolutionary Algorithm for Selective RNA Therapy of Trinucleotide Repeat Disorders" for more info
    
    Parameters:
        sensor_len (int)
        wt_len (int): Number of nucleotides of the wild-type repeats segment
        stride_size (int): Number of nucleotides to shift the WT repeats when moving from one WT hydrization state to the other
    
    Returns
        truth_vector (Boolean list): For each state, stores a 0 if we don't want the sRz to cleave 1 if we want it to cleave
    '''

    truth_vector = list()
    num_strides = int((sensor_len - wt_len) / stride_size) + 1

    # No binding. We don't want the ribozyme to cleave before the repeats have a chance to bind to the sensor
    truth_vector.append(0)

    # WT binding. We don't want to cleave the WT
    # Since the WT is shorter than the sensor, there are many different ways in which it can bind to the sensor
    # We consider one configuration (state) for each stride
    for s in range(num_strides):
        truth_vector.append(0)

    # Mutant binding. We want to cleave the mutant
    truth_vector.append(1)

    return truth_vector

def calculate_stride_size(num_repeats_wt, num_repeats_mutant, max_num_strides, repeat_length=3):
    '''
    The sensor is the same size as the mutant repeats segment, meaning that the wild-type repeats segment is shorter than the sensor.
    Therefore, the WT repeats can bind to the sensor at multiple possible sites.
    These sites are shifted by repeat_length nucleotides from each other.
    We refer to each of these possible shifts as a stride.
    The stride size is equal to repeat_length.
    The total number of strides is equal to (mutant length - wild length) / stride_size.
    For example, if the WT len is 30 and the mutant length is 60, there are (60-30)/3 = 10 possible strides
    However, considering each possible stride slows down fitness evaluation, since each stride is a folding state
    Therefore, we limit the maximum number of strides to max_num_strides.
    But, we want to make sure the strides are still spaced apart.
    This function calculates a new stride_size (which is an integer multiple of repeat_length) so that the max_num_strides is not exceeded
    '''

    return ceil( (num_repeats_mutant - num_repeats_wt) / float(max_num_strides) )*repeat_length

def generate_folding_tasks(sensor_mismatches, strands, name_to_bases, wt_len, stride_size):
    '''
    Generates a list of folding tasks to be executed by the folding software (in our case Vienna RNA)

    Parameters:
        sensor_mismatches (list of ints): 0-indexed positions of the mismatches w.r.t. start of the sensor segment
        strands (list of strings): strand 0 is the sRz, strand 1 is the (possibly extended) RzBS
        name_to_bases (dict): map from segment names to 1-indexed positions of each segment's nucletoides w.r.t. to the 5' end of the first strand
        wt_len (int)
        stride_size (int)
    
    Returns:
        tasks: A list containing one entry of the form (sequence, folding constraint) per hybridization state. 
               sequence is always the same, but folding constraint varies from state to state since it is used to simulate repeats binding to the sensor
    
    '''

    tasks = list()

    sensor_base_ids = name_to_bases['SENSOR'] # These are 1-indexed and right inclusive
    sensor_lims = (sensor_base_ids[0], sensor_base_ids[-1] + 1) # These are 1-indexed and right exclusive . ##?????????????????????????????????????????????????????????? DOES THIS MESS UP THE MISMATCH POSITIONS????????
    sensor_len = len(sensor_base_ids)
    strand_0_len = len(strands[0])
    strand_1_len = len(strands[1])

    num_strides = int((sensor_len - wt_len) / stride_size) + 1

    sequence = strands[0] + "&" +  strands[1]
    amp = '&'

    sensor_mismatches_absolute = [o + sensor_lims[0] for o in sensor_mismatches]

    # No binding. We don't want the ribozyme to cleave before the repeats have a chance to bind to the sensor
    fc = generate_folding_constraint_string(strand_0_len, [ ]) + amp + generate_folding_constraint_string(strand_1_len, [ ])
    tasks.append((sequence, fc))

    # WT binding. We don't want to cleave the WT
    # Since the WT is shorter than the sensor, there are many different ways in which it can bind to the sensor
    # We consider one configuration (state) for each stride
    folding_constraints = list()
    lims = [None, None]
    for s in range(num_strides):
        lims[0] = sensor_lims[0] + s*stride_size
        lims[1] = lims[0] + wt_len
        sensor_constraint = [i for i in range(lims[0], lims[1]) if i not in sensor_mismatches_absolute] 
        fc = generate_folding_constraint_string(strand_0_len, sensor_constraint ) + amp +  generate_folding_constraint_string(strand_1_len, [ ])
        tasks.append((sequence, fc))

    # Mutant binding. We want to cleave the mutant
    sensor_constraint = [i for i in range(sensor_lims[0], sensor_lims[1]) if i not in sensor_mismatches_absolute] 
    fc = generate_folding_constraint_string(strand_0_len, sensor_constraint ) + amp + generate_folding_constraint_string(strand_1_len, [ ])
    tasks.append((sequence, fc))

    return tasks
