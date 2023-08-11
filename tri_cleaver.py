from tri_cleaver_ea import TriCleaver

import argparse
import os

def read_file(filename):
    with open(filename, 'r') as file:
        return file.read()

def parse_args():
    parser = argparse.ArgumentParser(description='TriCleaver')

    # Required arguments
    parser.add_argument('--wt_file', type=str, required=True, help='The name of the .txt file containing the sequence of the wild-type strand of the target trinucleotide repeat disorder (TRD). The repeat region must be marked with a "-" on each side. To target one of the diseases studied in the accompaning paper, use "MDJ_wt.txt", "OPMD_wt.txt" or "HUNT_wt.txt". ')
    parser.add_argument('--mutant_file', type=str, required=True, help='The name of the .txt file containing the sequence of the mutant strand of the target trinucleotide repeat disorder (TRD). The repeat region must be marked with a "-" on each side. Regions upstream and downstream the repeat region must be identicial for the both wild-type and mutant strand. To target one of the diseases studied in the accompaning paper, use "MDJ_mutant.txt", "OPMD_mutant.txt" or "HUNT_mutant.txt". ')

    # Optional arguments
    parser.add_argument('--folding_temperature', type=float, default=37.0, help='Folding temperature')



    return parser.parse_args()

def validate_args(args):
    filenames = [args.wt_file, args.mutant_file]

    for i, filename in enumerate(filenames, 1):
        # Check if the provided filename has a .txt extension. If not, add it.
        if not filename.endswith(".txt"):
            filenames[i-1] += ".txt"

        # Check if the file exists
        if not os.path.exists(filenames[i-1]):
            raise ValueError(f"File '{filenames[i-1]}' not found!")

    args.wt_file, args.mutant_file = filenames
    return args

def count_triplets(sequence):
    offsets = [0, 1, 2]
    global_most_common_triplet = None
    global_max_count = 0
    for offset in offsets:
        triplets = [sequence[i:i+3] for i in range(offset, len(sequence), 3)]
        triplets = [t for t in triplets if len(t) == 3]
        triplet_counts = {}
        for triplet in triplets:
            if triplet not in triplet_counts:
                triplet_counts[triplet] = 0
            triplet_counts[triplet] += 1
        most_common_triplet = max(triplet_counts, key=triplet_counts.get)
        if triplet_counts[most_common_triplet] > global_max_count:
            global_most_common_triplet = most_common_triplet
            global_max_count = triplet_counts[most_common_triplet]
    
    return (global_most_common_triplet, global_max_count)

def process_input_seqs(wt_seq, mutant_seq, s1_len, s3_len):
    valid_chars = set("AUGTCaugtc-")    
    if not wt_seq: 
        raise ValueError(f"Please provide the sequence of the wild-type strand")
    if not set(wt_seq).issubset(valid_chars):
        raise ValueError("Wild-type sequence contains invalid characters. Allowed characters are A, U, G, T, C, as well a dash (-) on each side of the repeats segment")
    if not mutant_seq: 
        raise ValueError(f"Please provide the sequence of the wild-type strand")
    if not set(mutant_seq).issubset(valid_chars):
        raise ValueError("Mutant sequence contains invalid characters. Allowed characters are A, U, G, T, C, as well a dash (-) on each side of the repeats segment")
 
    #Extract the upstream, downstream, and repeat stretches
    if wt_seq.count('-') != 2 or mutant_seq.count('-') != 2:
        raise ValueError("Input sequences should each have exactly two dashes.")
    upstream_wt, repeats_wt, downstream_wt = wt_seq.split('-')
    upstream_mutant, repeats_mutant, downstream_mutant = mutant_seq.split('-')
    if downstream_wt != downstream_mutant:
        raise ValueError("The sequence of the region downstream of the repeats must be same in both the wild-type and mutant")
    downstream_seq = downstream_wt
    if upstream_wt != upstream_mutant:
        raise ValueError("The sequence of the region upstrea of the repeats must be same in both the wild-type and mutant")
    upstream_seq = upstream_wt
    #Determine the # of wt and mutant repeats
    triplet_wt, num_repeats_wt = count_triplets(repeats_wt)
    triplet_mutant, num_repeats_mutant = count_triplets(repeats_mutant)
    if triplet_wt != triplet_mutant:
        raise ValueError("The type of repeat (e.g. GCG, CAG) must be the same between the wild-type and mutant.")
    if num_repeats_wt >= num_repeats_mutant:
        raise ValueError("The number of mutant repeats must exceed the number of wild-type repeats.")
    if num_repeats_wt < 5:
        raise ValueError("The wild type must have at least 5 repeats.")
    repeat_unit = triplet_wt

    #In some cases (e.g. Hunts, MJD) there are a lot of mutant repeats relative to the wild-type. In one sense this is good: the higher the ratio between them, the easier it is for the OBS to differentiate between them. However, at some point there is overkill. The OBS will be very long which will grind folding to halt without breaking any more increase sensitivity. Therefore, we cap the OBS length.
    max_mutant_wt_difference = 20
    mutant_wt_ratio     = num_repeats_mutant / num_repeats_wt
    if num_repeats_mutant - num_repeats_wt > max_mutant_wt_difference:
        num_repeats_mutant = num_repeats_wt + max_mutant_wt_difference

    #Count the cut sites
    sbs_len = s1_len + 1 + s3_len
    cut_site_counts = dict()
    num_cut_sites = 0
    for i in range( len(downstream_seq) - sbs_len ):
        cut_site = downstream_seq[i + s3_len - 2:i + s3_len + 1]
        if cut_site[1] == 'U' and cut_site[2] == "G": #NUH
            num_cut_sites += 1
            if cut_site not in cut_site_counts.keys():
                cut_site_counts[cut_site] = 0
            else:
                cut_site_counts[cut_site] += 1
    if num_cut_sites == 0:
        raise ValueError("No potential cut sites detected. Unable to process request.")
    allowed_cut_site_types = ['GUC', 'AUC', 'UUC', 'AUA', 'GUA', 'CUC', 'UUA', 'AUU', 'GUU', 'CUA', 'CUU', 'UUU']


    return downstream_seq, num_repeats_wt, num_repeats_mutant, repeat_unit, allowed_cut_site_types


if __name__ == '__main__':
    #Harcoded parameters
    viability_params    = [-1.00, 0.90, 1,  0.20, 49]
    num_inds            = 300 #600
    num_gens            = 200
    mutation_rate       = 6
    mutation_probs      = [0.85, 0.15]
    selected_fitnesses  = ['noveltyCut', 'ONCut', 'OFFCut', 'mismatchCut']
    spacing         = 10
    spacing_length  = 2
    max_num_strides = 10
    s1_len, s3_len = 7,7
    sbs_len = s1_len + 1 + s3_len

    #Input parameters
    args = parse_args()

    try:
        args       = validate_args(args)        
        wt_seq     = read_file(args.wt_file)
        wt_seq = wt_seq.replace(' ','').replace('\n','').replace('a','A').replace('g','G').replace('c','C').replace('t','U').replace('T','U')
        mutant_seq = read_file(args.mutant_file)
        mutant_seq = mutant_seq.replace(' ','').replace('\n','').replace('a','A').replace('g','G').replace('c','C').replace('t','U').replace('T','U')
        downstream_seq, num_repeats_wt, num_repeats_mutant, repeat_unit, allowed_cut_site_types = process_input_seqs(wt_seq, mutant_seq, s1_len, s3_len)
    except ValueError as e:
        print(f"Error: {e}")
        exit(1)

    template_params = (spacing, spacing_length, s1_len, s3_len, sbs_len, max_num_strides, allowed_cut_site_types)

    bool_fast_folding, dot_dir = 1, 1 #dot_dir is vestigial
    job = TriCleaver(dot_dir, viability_params, bool_fast_folding, args.folding_temperature, downstream_seq, num_repeats_wt, num_repeats_mutant, repeat_unit, template_params)

    print("Generating selective ribozymes")
    job.run(num_inds=num_inds, num_gens=num_gens, mutation_rate=mutation_rate, mutation_probs=mutation_probs, selected_fitnesses=selected_fitnesses, parent_selection_method=0, parent_selection_params=[], survivor_selection_method=0, crossover_prob=0,  bool_store=1)
