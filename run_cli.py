import argparse
import os
import sys
import pickle

from ribogate.tc.substrate_processing import count_tandem_repeats, extract_every_rzbs
from ribogate.tc.tri_cleaver import TriCleaver
import config
from ribogate.validation import validate_wt_mutant_seqs, validate_repeats, validate_potential_cut_sites, validate_selected_cut_site_motifs, validate_environmental_parameters, validate_temperature
from ribogate.common.reporting import export_spreadsheet

def get_job_status(_):
    '''
    Only reason this function here is for consistency with server code that allows job interruptions
    '''
    return 'Running'

def save_checkpoint(self, _, checkpoint_rate):
    checkpoint_dir = "checkpoints"
    if (self.current_generation + 1) % checkpoint_rate == 0:

        if not os.path.exists(checkpoint_dir):
            os.makedirs(checkpoint_dir)

        checkpoint_file = os.path.join(checkpoint_dir, f'checkpoint_gen_{self.current_generation}.pkl')

        # Delete any existing checkpoint files
        for f in os.listdir(checkpoint_dir):
            if f.endswith('.pkl'):
                os.remove(os.path.join(checkpoint_dir, f))

        with open(checkpoint_file, 'wb') as f:
            pickle.dump(self.population, f)

def delete_all_checkpoints(checkpoint_dir = 'checkpoints'):   
    for f in os.listdir(checkpoint_dir):
        if f.startswith('checkpoint_gen_') and f.endswith('.pkl'):
            os.remove(os.path.join(checkpoint_dir, f))    

def load_checkpoint(checkpoint_dir = 'checkpoints'):   
    if not os.path.exists(checkpoint_dir):
        return None, -1

    checkpoints = [f for f in os.listdir(checkpoint_dir) if f.startswith('checkpoint_gen_')]
    if not checkpoints:
        return None, -1

    # Find the latest checkpoint file
    checkpoint_files = [f for f in os.listdir(checkpoint_dir) if f.endswith('.pkl')]
    if not checkpoint_files:
        raise FileNotFoundError('No checkpoint files found')

    latest_checkpoint = max(checkpoint_files, key=lambda f: int(f.split('_')[2].split('.')[0]))

    with open(os.path.join(checkpoint_dir, latest_checkpoint), 'rb') as f:
        population = pickle.load(f)
        # The population was saved at the end of this generation
        generation = int(latest_checkpoint.split('_')[-1].split('.')[0]) 

    return population, generation

def check_existing_checkpoint_and_prompt(checkpoint_dir='checkpoints'):
    # Check if the checkpoint directory exists and contains checkpoint files
    if os.path.exists(checkpoint_dir) and any(f.startswith('checkpoint_gen_') and f.endswith('.pkl') for f in os.listdir(checkpoint_dir)):
        print('A checkpoint already exists. Do you want to:')
        print('  1. Continue from the saved checkpoint')
        print('  2. Start a new run')
        choice = input('Enter your choice (1 or 2): ').strip()

        if choice == '1':
            population, generation = load_checkpoint(checkpoint_dir)
            return population, generation
        elif choice == '2':
            delete_all_checkpoints(checkpoint_dir)             
            return None, -1  # Start from generation 0
        else:
            print('Invalid choice. Exiting...')
            sys.exit(1)
    else:
        return None, -1  # Start from generation 0


def read_file(filename):
    with open(filename, 'r') as file:
        return file.read()

def parse_args():
    parser = argparse.ArgumentParser(description='TriCleaver')

    # Required arguments
    parser.add_argument('--wt_file', type=str, required=True, help='The name of the .txt file containing the sequence of the wild-type strand of the target trinucleotide repeat disorder (TRD). The repeat region must be marked with a "-" on each side. To target one of the diseases studied in the accompaning paper, use "MJD_wt.txt", "PAPBN1_wt.txt" or "HTT_wt.txt". ')
    parser.add_argument('--mutant_file', type=str, required=True, help='The name of the .txt file containing the sequence of the mutant strand of the target trinucleotide repeat disorder (TRD). The repeat region must be marked with a "-" on each side. Regions upstream and downstream the repeat region must be identicial for the both wild-type and mutant strand. To target one of the diseases studied in the accompaning paper, use "MJD_mutant.txt", "PABPN1_mutant.txt" or "HTT_mutant.txt". ')

    # Optional arguments
    parser.add_argument('--vienna-path', type=str, required=False, default='vienna_rna_binaries')
    parser.add_argument('--folding-temperature', type=float, default=37.0, help='Folding temperature')
    parser.add_argument('--selected-cut-site-motifs', type=str, default='GUC,AUC,UUC,AUA,GUA,CUC,UUA,AUU,GUU,CUA,CUU,UUU', help='A comma seperated list of allowed cut site motifs for the hammerhead ribozyme. These motifs must be of the form NUH. Example:CUA,CUG,CUA,GUC')
    parser.add_argument('--environment', type=str, choices=['in_vivo_alen', 'in_vivo_kharma', 'in_vivo_ribozyme_only','in_vitro_ribozyme_only','in_vitro_default_T7'], default='in_vivo_ribozyme_only', help='The environment in which the selective ribozymes (sRzs) will be tested. Currently this has no effect on the actual design of the sRzs. Instead, it determines how the sRzs will be packaged in preparation of DNA ordering')
    parser.add_argument('--num-processes', type=int, required=False, default=12) 

    #args = parser.parse_args()    
    args, unknown = parser.parse_known_args() # This makes it easier for debugging the CLI and web app from a single set of command arguments
    
    return args

def raise_error_from_dict(errors_dict):
    # Raise the first error encountered in the dict
    for error_type, errors in errors_dict.items():
        if len(errors) > 0:
            raise ValueError(errors[0])

def validate_filenames(wt_file, mutant_file):
    
    # Validate the filenames
    filenames = [wt_file, mutant_file]
    for i, filename in enumerate(filenames, 1):
        # Check if the provided filename has a .txt extension. If not, add it.
        if not filename.endswith(".txt"):
            filenames[i-1] += ".txt"

        # Check if the file exists
        if not os.path.exists(filenames[i-1]):
            raise ValueError(f"File '{filenames[i-1]}' not found!")
        
    return filenames


if __name__ == '__main__':
    args = parse_args()

    # Process wt and mutant sequences
    wt_file, mutant_file = validate_filenames(args.wt_file, args.mutant_file)
    wt_seq = read_file(wt_file)
    mutant_seq = read_file(mutant_file)    
    repeat_unit_len = 3 # For now only allow Trinucleotide repeats
    data, errors = validate_wt_mutant_seqs(wt_seq, mutant_seq)
    if errors != None:
        raise_error_from_dict(errors)

    wt_seq = data['wt_seq']
    mutant_seq = data['mutant_seq']
    downstream_seq = data['downstream_seq']
    upstream_seq = data['upstream_seq']
    repeats_wt = data['repeats_wt']
    repeats_mutant = data['repeats_mutant']
    
    # Process repeats
    data, errors = validate_repeats(repeats_wt, repeats_mutant, repeat_unit_len)
    if errors != None:
        raise_error_from_dict(errors)
            
    num_repeats_wt = data['num_repeats_wt']
    num_repeats_mutant = data['num_repeats_mutant']
    repeat_unit = data['repeat_unit']
    
    # Process ribozyme type
    ribozyme_type = 'minimal_hammerhead'
    
    # Process potential cut sites
    data, errors = validate_potential_cut_sites(ribozyme_type, downstream_seq)
    if errors != None:
        raise_error_from_dict(errors)
            
    rzbs_positions = data['rzbs_positions']
    rzbs_cut_site_motifs = data['rzbs_cut_site_motifs']     
    num_cut_site_motifs = data['num_cut_site_motifs']
    
    # Process temperature            
    data, errors = validate_temperature(args.folding_temperature)
    if errors != None:
        raise_error_from_dict(errors)
            
    temperature = data['temperature']            

    # Process environmental parameters
    environment = args.environment                               
  
    # Process selected cut site motifs
    data, errors = validate_selected_cut_site_motifs(rzbs_positions, rzbs_cut_site_motifs, args.selected_cut_site_motifs.split(',') )
    if errors != None:
        raise_error_from_dict(errors)
            
    selected_cut_site_motifs = data['selected_cut_site_motifs']
    allowed_rzbs_positions = data['allowed_rzbs_positions']
    
    # Set vienna path
    config.folding_package_paths = {'vienna_binaries_max_loop_30':args.vienna_path, 'vienna_binaries_max_loop_300':args.vienna_path, 'vienna_api':None}
    
    # Set BLAST parameters. For now always use False since we don't require stand-alone users to have a local blast database
    config.should_perform_blast = False
    config.blast_db_path = None
    
    # Set number of processes
    config.num_processes = args.num_processes

    # Always set this to this False.
    config.should_report_unhybridized_srz = False
   
    # Create TriCleaver EA instance    
    job = TriCleaver(ribozyme_type, config.template_params_tc, args.folding_temperature, downstream_seq, num_repeats_wt, num_repeats_mutant, repeat_unit, allowed_rzbs_positions)
    
    # Set the general EA parameters
    job_id = None # This is only needed when launched from the web server
    job.set_ea_params(config.ea_params_tc, save_checkpoint, get_job_status, job_id)

    # Check if there were any saved checkpoints
    initial_population, last_completed_generation = check_existing_checkpoint_and_prompt()

    # Run the EA
    print("Generating selective ribozymes")
    job.run(initial_population, last_completed_generation)

    # Once the job has completed, delete all checkpoints to avoid confusing users when they want to create the next new job
    delete_all_checkpoints()    

    export_spreadsheet(job.population, {'ENVIRONMENT':environment}, 'tc', 'tricleaver_results.xlsx') 