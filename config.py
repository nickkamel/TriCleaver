'''
Defines parameter values that are used in both the web server and CLI versions
Web app specific parameter values are stored in the app.config dictionary which is populated in web_app.py and runserver.py
SRz and ribogate template parameters are currently defined in ribozyme_templates.py
Contains parameters for both TriCleaver (TC) and TruthSeqEr (TS)
The TS parameters are maintained in the TC CLI code, but serve no purpose
'''

import importlib.util

# We can either treat the core as one single motif, or 3 motifs (C1, C2, C3) when evaluating ribozyme activity
# Using a single motif is cleaner and I don't think there is much of a difference between the two cases
should_treat_core_as_single_motif = True

# This is monospace meaning that the sequence, constraint, and dot bracket are perfectly aligned
spreadsheet_font_name = 'Consolas' 

# Limits for the temperature the user can specify
min_temp = 20
max_temp = 80

####### TC ribozyme template parameters ##########
template_params_tc = dict()

# Length of the linker segments connecting the sensor to the stem 2
template_params_tc['LINKER_LEN'] = 7

# Stem lengths for the default ribozyme model (minimal hammerhead)
template_params_tc['STEM_LENS_MINIMAL'] = {'S1':7, 'S2':7, 'S3':7}

# Stem lengths for the 'extended' hammerhead
# Note that this hammerhead isn't really the extended version since we disrupt stem 2 to add the extension region
# Therefore, it is currently used for research purposes only
# However, it does provide a roadmap for adding other ribozymes
# Stem 1 of the extended hammerhead contains a bulge loop, so I use the non-conventional naming stem0 and stem1 to designate the two sides of this loop
template_params_tc['STEM_LENS_EXTENDED'] = {'S0':4, 'S1':5, 'S2':7, 'S3':7}

# To prevent long stretches of double-stranded RNA when the repeats binding to the sensor, we add NUM_CONTIGUOUS_MISMATCHES mismatches between the sensor-repeats every 'MISMATCH_SPACING' nucleotides
template_params_tc['NUM_CONTIGUOUS_MISMATCHES'] = 2
template_params_tc['MISMATCH_SPACING'] = 9

# Determines the maximum number of sRz + Mutant folding configurations 
# This has a major effect on performance since the greater the nubmer of strides, the greater the number of hydridization states that are folded
template_params_tc['MAX_NUM_STRIDES'] = 10 

# During the EA, extend the RzBS by a short amount. After the EA, (optionally) extend it by a longer amount.
template_params_tc['RZBS_EXTENSION_LENS'] = {'LIVE': 10, 'POST':10} 
#template_params_tc['RZBS_EXTENSION_LENS'] = {'LIVE': 10, 'POST':50} # Longer Post extension but reduces the nubmer of viable individuals

# Cap if the difference in repeats between the mutant and WT is very large, we can cap the number of mutant repeats to reduce the size of sensor and accelerate folding
template_params_tc['MAX_MUTANT_WT_DIFFERENCE'] = 20

# This is used for coarse-graining the secondary structures
# If True, we consider the individual segments of the ribozyme (e.g. S1A, C1, S2A, etc.)
# If False, we consider treat the ribozyme as a single segment
# I recommend always setting this to True
template_params_tc['ARE_SEGMENTS_LOW_LEVEL'] = True

# Determines how much importance we give to two individuals having different RzBSs when measuring their phenotypic distance
# Set to 0 to only consider the structural difference between two individuals
# Set to 1 to scale the structural difference by 2 when the two individuals have distinct RzBS. 
# Setting it to 1 strikes a good balance between emphasizing structural diversity and RzBS diversity
template_params_tc['RZBS_DIVERSITY_WEIGHT'] = 1
##################################################

# Check to see whether user has the Vienna RNA python bindings installed
# They can be installed with Anaconda, but I have never been able to get them to work with just pip
# They must required some non-default build tools 
# Therefore, I don't want to require the paper reviewers to install the Vienna bindings
# Also, the bindings seem to have a hardcoded 30 MaxLoop
# So if they are installed, we use those in lieu of the default vienna binaries that are compile with 30 MaxLoop
if importlib.util.find_spec('RNA') is not None:
    vienna_fast_package_name = 'vienna_api'
else:
    vienna_fast_package_name = 'vienna_binaries_max_loop_30'

# TC folding parameters
# Use the smaller loop size (30) during the EA (LIVE) and the larger loop size (300) during post-processing (POST)    
folding_packages_tc = dict()
folding_packages_tc['LIVE'] = [vienna_fast_package_name]
folding_packages_tc['POST'] = ['vienna_binaries_max_loop_300']

# TC EA params
viability_params_tc = {'MIN_VIA':0, 'MAX_VIA':0.95, 'HAS_BREAKPOINT':True, 'BREAKPOINT_VALUE':0.60, 'BREAKPOINT_GEN':49} # New params: new viability val = (old val + 1)/2 #Old params: [-1.00, 0.90, 1,  0.20, 49]   
ea_params_tc = {'NUM_INDS':300, 'NUM_GENERATIONS':200, 'MUTATION_RATE':6, 'MUTATION_PROBABILITIES':[0.85, 0.15], 'CROSSOVER_PROBABILITY':0, 'PARENT_SELECTION_METHOD': 'UNIFORM', 'PARENT_SELECTION_PARAMS':None,  'SURVIVOR_SELECTION_METHOD':'GENITOR', 'NOVELTY_NEIGHBORHOOD_SIZE':30, 'SELECTED_FITNESSES':['novelty', 'ON', 'OFF', 'mismatch'], 'IS_VIABILITY_NULLIFICATION_ENABLED':True, 'VIABILITY_PARAMS':viability_params_tc} # use 'NUM_INDS':8, 'NUM_GENERATIONS':4 for quick debugging

# Set the criteria for results of the final population to be reported to the user
# Multiply by 0.88 to give a bit of breathing room since we fold the sRz with largest loop size only in the post-processing step, which might result in some discrepencies
# I put 0.88 because all of the designs generated for the Nature Comms Revision had a viability above 0.88
# When using the extended hammerhead, we might want to lower to floor for reported viability.
reportability_criteria_tc = [{'name':'viability', 'floor':0.88, 'ceil':float('inf')}]

# TS folding parameters
# For some 3-input functions like f-105, fictious mechanisms are discovered if we use the smaller loop size during the run
# Therefore, the larger loop size should be used during the LIVE run. 
# But this slows down debugging, so for now we just use it for POST processing
folding_packages_ts = dict()
folding_packages_ts['LIVE'] = ['vienna_binaries_max_loop_300']
#folding_packages_ts['LIVE'] = [vienna_fast_package_name]
folding_packages_ts['POST'] = ['vienna_binaries_max_loop_300']

# This determines whether we need to import the RNA module inside folding.py
# I'm no longer using this. I am just import RNA directly from the folding function. This is less ideal but minimizes inadvertent edge cases
#is_vienna_api_used = 'vienna_api' in folding_packages_tc['LIVE'] or 'vienna_api' in folding_packages_tc['POST'] or 'vienna_api' in folding_packages_ts['LIVE'] or 'vienna_api' in folding_packages_ts['POST']

# TS EA params
viability_params_ts = {'MIN_VIA':0, 'MAX_VIA':0.95, 'HAS_BREAKPOINT':True, 'BREAKPOINT_VALUE':0.475, 'BREAKPOINT_GEN':49} # New params: new viability val = (old val + 1)/2 #Old params: [-1.00, 0.90, 1, -0.05, 49]
ea_params_ts = {'NUM_INDS':300, 'NUM_GENERATIONS':200, 'MUTATION_RATE':4, 'MUTATION_PROBABILITIES':[1, 0], 'CROSSOVER_PROBABILITY':0, 'PARENT_SELECTION_METHOD': 'UNIFORM', 'PARENT_SELECTION_PARAMS':None,  'SURVIVOR_SELECTION_METHOD':'GENITOR', 'NOVELTY_NEIGHBORHOOD_SIZE':30, 'SELECTED_FITNESSES':['novelty', 'ON', 'OFF'], 'IS_VIABILITY_NULLIFICATION_ENABLED':True, 'VIABILITY_PARAMS':viability_params_ts}

# We can later add a criterion for verification as well if need be
reportability_criteria_ts = [{'name':'viability', 'floor':0.88, 'ceil':float('inf')}] 

####### TS ribozyme template parameters ##########
template_params_ts = dict()
template_params_ts['STEM_LENS'] = {'S1':7, 'S2':8, 'S3':7}

# If we need to disrupt the ribozyme when all OBSs are occupied by inputs, we require an additional segment (the logic linker)
template_params_ts['LOGIC_LINKER_LEN'] = 14

# We may want to provide some flexibility between OBSs
template_params_ts['HAS_FLEXIBLE_LINKERS'] = True
template_params_ts['FLEXIBLE_LINKER_LEN'] = 3

# Set this to true if we want to transcribe the ribogates in vitro
template_params_ts['HAS_STARTING_GG'] = True

# See explanation in the TC template params section
# I used False to generate all the designs for mechanism extraction
template_params_ts['ARE_SEGMENTS_LOW_LEVEL'] = True 
################################################

# Colors of the selective ribozymes (designed by TC) and ribogate (designed by TS) segments
s1_color = '#e6194b'
s2_color = '#ffe119' 
s3_color = '#000000'
core_color = '#3cb44b'
linker_color = '#4363d8'
sensor_color = '#f58231'
rzbs_extension_color = '#888888'
pseudoknot_color = '#888888' #for extended harmmerhead
obs1_color = '#f58231'
obs2_color = '#46f0f0'
obs3_color = '#BCF60C'
s3_hairpin_color = '#888888'
repeats_color = '#911eb4' # Repeats are never actual drawn
seg_colors = {'S1A': s1_color, 'S1B': s1_color, 'S2A': s2_color, 'S2B': s2_color, 'S3A': s3_color, 'S3B': s3_color, 'C1': core_color,  'C2': core_color, 'C3': core_color, 
                                'L0': linker_color, 'L1': linker_color, 'SENSOR': sensor_color, 'UPSTREAM_RZBS': rzbs_extension_color, 'DOWNSTREAM_RZBS': rzbs_extension_color, 'S0A': s1_color, 'S0B': s1_color, 'S1KNOT': pseudoknot_color, 'S2KNOT': pseudoknot_color, 'OBS1':obs1_color, 'OBS2':obs2_color, 'OBS3':obs3_color, 'L2': linker_color, 'L3': linker_color, 'S3H':s3_hairpin_color, 'startGG':s1_color}
verbose_seg_names = {'S1A': '1st half of stem 1', 'S1B': '2nd half of stem 1', 'S2A': '1st half of stem 2', 'S2B': '2nd half of stem 2', 'S3A': '1st half of stem 3', 'S3B': '2nd half of stem 3', 'C1': '1st part of core', 'C2': '2nd part of core', 'C3': '3rd part of core', 'L0': '1st linker', 'L1': '2nd linker', 'SENSOR': 'sensor', 'UPSTREAM_RZBS': 'region upstream of RzBS', 'DOWNSTREAM_RZBS': 'region downstream of RzBS', 'S0A': "1st half of stem 1", 'S0B': "2nd half of stem 1", 'S1KNOT': "1st part of pseudoknot", 'S2KNOT': "2nd part of pseudoknot", 'OBS1':"Oligoonucleotide binding site for input 1",  'OBS2':"Oligoonucleotide binding site for input 2",  'OBS3':"Oligoonucleotide binding site for input 3", 'L2':'3rd linker', 'L3':"4th linker", 'S3H': "Stem 3 hairpin", 'startGG':"5' GG"}  

# Constructs
# The selective ribozymes designed by TC may be incorporated constructs with other segments for DNA ordering
construct_segment_names = dict()
construct_segment_names['in_vivo_alen'] = ["partial tRNA val human", "linker", "selective ribozyme", "restriction enzyme cut site"] # Use this when splicing in the sRz into the Nawrot plasmid
construct_segment_names['in_vivo_kharma'] = ["tRNA val human", "selective ribozyme", "CTE"] # Use this for the new plasmid
construct_segment_names['in_vivo_ribozyme_only'] = ["selective ribozyme"] # Only generate the pure selective ribozyme sequence
construct_segment_names['in_vitro_ribozyme_only'] = ["selective ribozyme"] # Only generate the pure selective ribozyme sequence
construct_segment_names['in_vitro_default_promoter'] = ["T7", "selective ribozyme"] # Add a T7 promoter for in vitro transcription
construct_segment_names['in_vitro_custom_promoter'] = ["T7", "selective ribozyme"] # Add a custom promoter for in vitro transcription

construct_segments = dict()
construct_segments['partial tRNA val human'] = {'color': 'red', 'description':"Portion of human tRNA to faciliate the export of the selective ribozyme from the nucleus", 'sequence': "TTCGAAACCGGGCACTACAAAAACCAAC"}
construct_segments['tRNA val human'] = {'color':'red', 'description':"Human tRNA to faciliate the export of the selective ribozyme from the nucleus", 'sequence': "GTTTCCGTAGTGTAGTGGTTATCACGTTCGCCTAACACGCGAAAGGTCCCCGGTTCGAAACCGGGCGGAA"}
construct_segments['linker'] = {'color': 'black', 'description': "Linker sequence", 'sequence':"TTT"}
construct_segments['selective ribozyme'] = {'color':'black', 'description': "The selective ribozyme designed by TriCleaver", 'sequence':None} # We don't know the sequence of the sRz a priori (it is designed by TriCleaver)
construct_segments['restriction enzyme cut site'] = {'color':'blue', 'description':"Cut site to help ligate the designed selective ribozyme into the Nawrot plasmid", 'sequence':"GGTAC"}
construct_segments['CTE'] = {'color':'green', 'description':"Constituative transport element which serves as a binding site for a helicase that unwinds the substrate to faciliate ribozyme binding", 'sequence':"AGACCACCTCCCCTGCGAGCTAAGCTGGACAGCCAATGACGGGTAAGAGAGTGACATTGTTCACTAACCTAAGACAGGAGGGCCGTCAGAGCTACTGCCTAATCCAAAGACGGGTAAAAGTGATAAAAATGTATCACTCCAACCTAAGACAGGCGCAGCTTCCGAGGGATTTG"}
construct_segments['promoter'] = {'color':'red', 'description':"T7 promoter", 'sequence':"GAATTTAATACGACTCACTATA" } # An appropriate number of Gs will be dynamically added once the ribozyme sequence is known