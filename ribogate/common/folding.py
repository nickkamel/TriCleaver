'''
Interfaces with a customized version of Vienna RNA.
This code could be modified to interface with other folders, so look as they output 1) a base-pairing probability matrix (BPPM) and 2) a MFE dot-bracket structure
'''

import subprocess
from collections import defaultdict

def fold(tasks, package_names, folding_temperature, num_strands, package_paths):
    all_package_results = list()
    num_packages = len(package_names)
    for i, package_name in enumerate(package_names):        
        if package_name == 'vienna_binaries_max_loop_30':
           package_results = fold_using_vienna_binaries(tasks, False, folding_temperature, num_strands, package_paths[package_name]) 
        elif package_name == 'vienna_binaries_max_loop_300':
           package_results = fold_using_vienna_binaries(tasks, True, folding_temperature, num_strands, package_paths[package_name])  
        elif package_name == 'vienna_api':
           package_results = fold_using_vienna_api(tasks, folding_temperature, num_strands) 
        else:
            raise NotImplementedError(f"There is no folding package named: ({package_name})")
        all_package_results.append(package_results)

    if num_packages == 1:
        # If only one folding package was used, just return its results
        return all_package_results[0]
    else:
        '''
        TO DO: If multiple folding packages are used, the BPPM of each package need to be combined somehow into a single BPPM
        Same for MFE structures
        '''
        pass


def fold_using_vienna_api(tasks, folding_temperature, num_strands):
    '''
    Folds a batch of (sequence, structural constraint) pairs using the Vienna RNA Python bindings

    Returns:
        processed results (list): One entry per (sequence, structural constraint) pair that was foldeded
            mfe_db (string) Minimum free energy structure in dot-bracket notation
            bppm (dict) Base-pairing probabiliy matrix encoding the probability of any two nucleotides being paired
            mfe_ad: OBSOLETE
            mfe_energy (float): The free energy of the MFE structure
            ens_div (float): The diversity of the ensemble of secondary structures 
            
    Cleaner than having to call the vienna binaries, but the value of MaxLoop seems to be hardcoded to 30 which can be issue in certain situations
    '''
    
    # I don't like importing this module directly in the function, but I want to be sure it won't accidently be imported for users that aren't using the API
    import RNA

    input_data = ''
    for task in tasks:
        seq = task[0]
        constraint = task[1]
        input_data = input_data + seq + "\n" + constraint + "\n"
    input_data = input_data[:-1]  # Get rid of the last newline

    # It seems like we can read RNA.MAXLOOP, but not set it (setting RNA.MAXLOOP = value, changes it, but it has no downstream effect)
    # Unless we found a workaround, it means we will either have to continue distributing the binaries or recompile and distribute the Python bindings 
    # print("Max loop", RNA.MAXLOOP)
  
    # Model details
    md = RNA.md()
    md.uniq_ML = 1
    md.temperature = folding_temperature

    processed_results = []

    for task in tasks:
        seq = task[0]
        constraint = task[1]

        # Create fold_compound with constraints
        fc = RNA.fold_compound(seq, md)

        # Split constraints if there are multiple strands 
        # Also keep track of the split location (this assumes that there is a maximum of two strands being co-folded)
        amp_index = constraint.find('&')
        constraints = constraint.split('&')

        # fc.pf() doesn't seem to accept the '&' character. So we remove it from the constraint string
        # TO DO: Confirm that this doesn't cause an off by one error for any constraints specified after the '&'
        fc.constraints_add(''.join(constraints), RNA.CONSTRAINT_DB_DEFAULT)
        
        # Apparently fc.mfe() works for both single and multiple strands
        (mfe_structure, mfe_energy) = fc.mfe()

        # The folding output omits the '&' from the predicted MFE structure
        if num_strands == 2:
            mfe_structure = mfe_structure[:amp_index] + '&' + mfe_structure[amp_index:]     

        # Calculate partition function and base-pairing probability matrix (BPPM)
        fc.pf()                
        probs = fc.bpp()
        
        
        # Store base pair probabilities
        # The BPPM generated through the API is triangular
        bppm = defaultdict(int)
        adj_dict = dict() # This is used to help us efficiently calculate the probability of base being unpaired       
        num_bases = len(task[0]) # This includes the '&' for co-folding but it doesn't break anything
        
        for i in range(num_bases):
            adj_dict[i + 1] = list()
        for i in range(1, num_bases):
            for j in range(i+1, num_bases):
                if probs[i][j] > 1e-5:
                    bppm[(i,j)] = probs[i][j]             
                    bppm[(j,i)] = probs[i][j]   
                    
                    # Used for help calculating unpaired probability
                    adj_dict[i].append(j)
                    adj_dict[j].append(i)

        # Calculate and store the probability of each base being unpaired
        for base_id1 in adj_dict.keys():
            paired_sum = 0
            for base_id2 in adj_dict[base_id1]:
                paired_sum += bppm[(base_id1, base_id2)]
            unpaired_prob = 1 - paired_sum

            bppm[(base_id1, None)]  = unpaired_prob
            bppm[(None, base_id1 )] = unpaired_prob
        

        mfe_ad = None # TO DO: Remove this everywhere

        if num_strands == 1:            
            ens_div = fc.ensemble_defect(mfe_structure)   
            processed_results.append((mfe_structure, bppm, mfe_ad, mfe_energy, ens_div))
        else:
            # TO DO: See if ensemble_defect now works for dimers
            processed_results.append((mfe_structure, bppm, mfe_ad, mfe_energy))

    return processed_results

def fold_using_vienna_binaries(tasks, should_fold_large_loops, folding_temperature, num_strands, vienna_path):    
    '''
    Folds a batch of (sequence, structural constraint) pairs using the Vienna RNA binaries
    
    Returns
        processed results (list): One entry per (sequence, structural constraint) pair that was foldeded
            mfe_db (string) Minimum free energy structure in dot-bracket notation
            bppm (dict) Base-pairing probabiliy matrix encoding the probability of any two nucleotides being paired
            mfe_ad: OBSOLETE
            mfe_energy (float): The free energy of the MFE structure
            ens_div (float): The diversity of the ensemble of secondary structures
    
    Calls a modified version of Vienna via subprocess and parses the results that would be outputted to the console.
    Normally, Vienna generates a dot.ps file storing the square root of the probabibility of two nucleotides being paired.
    Only pairs of nucleotides with a pairing probability exceeding a certain threshold are saved the file.
    However, writing all these dot.ps files to the disk and reading them becomes a massive bottleneck.
    Therefore, I modified the C source code of Vienna to output these probabilites to the console instead of the file.
    This was an ugly hack but it worked.
    I also added delimters such as "BPPM_START" and "BPPM_END" to the Vienna source to make parsing the output easier.

    I wanted to define package_paths as a global variable in app.config, 
    This seems to initially work, but after a while app.config['package_paths'] becomes undefined. 
    I think the cause is the multi-processing, but its difficult to debug.
    Therefore I pass package_paths as a parameter, which is uglier. 

    If 1 strand is selected, RNAFold is used.
    If 2 strands are selected, RNACofold is used.

    One issue with Vienna is that that by default loops are limited to 30 nucleotides.
    This causes the EA to design fictitious ribogates by exploiting this limitation.
    Therefore, I changed the MAXLOOP parameter to 300 in energy_const.h and recompiled the source code
    However, this version is a few times slower than the MAXLOOP = 30.
    I'm not sure if this is because largers loops are slower to fold or because I compiled the source in a suboptimal way.
    Therefore, I have do versions of RNAFold and RNACofold: the vanilla version and the '300' version.
    '''

    input_data = ''
    for task in tasks:
        seq = task[0]
        constraint = task[1]
        input_data = input_data + seq + "\n" + constraint + "\n"
    input_data = input_data[:-1] # Get rid of the last newline

    if should_fold_large_loops:
        str300   = '300'
    else:
        str300   = ''
    if num_strands == 1:
        app_name = "\RNAfold"
    else:
        app_name = "\RNACofold"

    p1 = subprocess.Popen([vienna_path + app_name + str300 + ".exe","-p","-C","--noPS","--bppmThreshold=1e-5", "--temp=" + str(folding_temperature)], stdin=subprocess.PIPE,stdout=subprocess.PIPE)   
    
    input_data = str.encode(input_data)
    results = p1.communicate(input=input_data)[0]
    #results_lines = results.split('\r\n')
    results_lines = results.decode().split('\r\n')
    results_lines = results_lines[:-1] # Get rid of the last newline

    results_starts = list()
    for line_idx, line in enumerate(results_lines):
        if line == "Secondary structure":
            results_starts.append(line_idx)
    results_starts.append(len(results_lines)) # Adds a virtual start id for what would be the next result

    num_results = len(results_starts) - 1

    processed_results = list()
    for i in range(num_results):
        result_lines = results_lines[results_starts[i]:results_starts[i+1]]
        #print(result_lines)

        mfe_db = result_lines[1]
        mfe_ad = None # TO DO: Remove this everywhere
        mfe_energy = float(result_lines[3])
        if num_strands == 1:
            ens_div = float(result_lines[-1])

        for line_idx, line in enumerate(result_lines):
            if line == "BPPM start":
                bppm_start_line_idx = line_idx + 1
            if line == "BPPM end":
                bppm_end_line_idx = line_idx - 1

        #print(bppm_start_line, bppm_end_line)
        bppm = defaultdict(int)
        adj_dict = dict()
        num_bases = len(tasks[i][0]) 
        for i in range(num_bases):
            adj_dict[i + 1] = list()
        for line_idx in range(bppm_start_line_idx, bppm_end_line_idx + 1):
            line = result_lines[line_idx]
            split = line.split('\t')
            base_id1 = int(split[0])
            base_id2 = int(split[1])
            prob = float(split[2]) ** 2 # The values given in the ps file are the sqrt of the probability
            bppm[(base_id1, base_id2)] = prob
            bppm[(base_id2, base_id1)] = prob
            adj_dict[base_id1].append(base_id2)
            adj_dict[base_id2].append(base_id1)
        for base_id1 in adj_dict.keys():
            paired_sum = 0
            for base_id2 in adj_dict[base_id1]:
                paired_sum += bppm[(base_id1, base_id2)]
            unpaired_prob = 1 - paired_sum

            bppm[(base_id1, None)]  = unpaired_prob
            bppm[(None, base_id1 )] = unpaired_prob

        if num_strands == 1:
            processed_results.append((mfe_db, bppm, mfe_ad, mfe_energy, ens_div))
        else:
            processed_results.append((mfe_db, bppm, mfe_ad, mfe_energy)) # RNA CoFold doesn't output an ensemble diversity value
    return processed_results

