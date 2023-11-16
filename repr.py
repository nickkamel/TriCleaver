from folding import generate_fold_constraint, build_map, allowed_lut
import numpy as np
from utility import arg_sort, inv_list
import copy

#Selective ribozyme class
class SRS():
    def __init__(self):
        self.sccs     = list()
        self.segments = list() #This is always sorted
        self.strands  = list()

        #For constraints, core score. 
        #input: name of a segment
        #output: list of base (global) ids contained in that segment
        self.name_to_bases        = dict() 
        #For coarse graining. 
        #input:   global id of base 
        ##output: pair of (global) id of segment containing that base and position of base within that segment
        self.base_to_seg          = dict()  

        #bcc: base connected component
        #scc: segment connected component
        #bcc_to_scc: inputs a scalar integer representing the global id of a bcc. Outputs a pair of integer representing the scc that contains it and its local id within that scc
        self.bcc_to_scc          = dict()
        self.mutation_prob       = list()
        self.geno                = ''
        self.folding_constraints = list()
        self.potential_fitnesses = dict()
        self.folded_strand_ids   = [0, 1] #This is used for building the color string. The ribozyme and the substrate, not the repeats. 

    def update_representation(self, downstream_seq, sensor_mismatches, repeats, misc_lens, bool_init): #!!!!! THIS NEEDS TO BE A METHOD OF INDIVIDUAL, NOT EA
        s1_len, s2_len, s3_len, linker_len, rzbs_len, rzbs_extension_len  = misc_lens
        pre_rzbs  = downstream_seq[self.rs-rzbs_extension_len:self.rs]
        rzbs      = downstream_seq[self.rs:self.rs + rzbs_len]
        rzbs_s3b  = rzbs[0:s3_len]
        rzbs_c3   = rzbs[s3_len:s3_len+1]
        rzbs_s1b  = rzbs[s3_len+1:rzbs_len]
        post_rzbs = downstream_seq[self.rs + rzbs_len:self.rs + rzbs_len + rzbs_extension_len]

        #If this is the first generation, no bccs have been created yet. The SCC initialization will fill sample them with bases / base-pairs from the corresponding set of valid nucleotide assigments. 
        #For any other generation, we keep track of the bccs at the end of the previous generation and supply them to the SCC initialization. If these old bccs violate any new constraints (which may happen if the RzBS is mutated), new bccs are sampled.
        num_sccs = 11
        if bool_init:
            bccs = [None for i in range(num_sccs)]
        else:
            bccs = [scc.bccs for scc in self.sccs]

        self.sccs     = [None]*num_sccs
        self.sccs[0]  = SCC(('S1A', 'S1B'), [('N',rzbs_s1b[::-1][i]) for i in range(s1_len)], s1_len, [0,1], [0,3], ['#e6194b', '#e6194b'], bccs[0])
        self.sccs[1]  = SCC(('C1',), list('CUGAUGA'), 7, [0], [1], ['#3cb44b'], bccs[1])
        self.sccs[2]  = SCC(('S2A', 'S2B'), [('G', 'C')] + [('N','N') for i in range((s2_len-1))], s2_len, [0,0], [2,6], ['#4363d8', '#4363d8'], bccs[2])
        self.sccs[3]  = SCC(('L0',), ['N' for i in range(linker_len)], linker_len, [0], [3], ['#ffe119'], bccs[3])
        self.sccs[4]  = SCC(('SENSOR','REPEATS'), [('N', base) for base in repeats[::-1]], len(repeats), [0,2], [4,0], ['#f58231', None], bccs[4], mismatches = copy.deepcopy(sensor_mismatches))        
        self.sccs[5]  = SCC(('L1',), ['N' for i in range(linker_len)], linker_len, [0], [5], ['#ffe119'], bccs[5])
        self.sccs[6]  = SCC(('C2',), list('GAA'), 3, [0], [7], ['#3cb44b'], bccs[6])
        self.sccs[7]  = SCC(('S3A', 'S3B'), [('A','U')] + [('N',rzbs_s3b[::-1][i]) for i in range(1, s3_len)], s3_len, [0,1], [8, 1], ['#911eb4', '#911eb4'], bccs[7])
        self.sccs[8]  = SCC(('PRE_RZBS',), list(pre_rzbs), len(pre_rzbs),[1], [0], ['#FFFFFF'], bccs[8])
        self.sccs[9]  = SCC(('C3',), list(rzbs_c3), 1, [1], [2], ['#3cb44b'], bccs[9])
        self.sccs[10] = SCC(('POST_RZBS',), list(post_rzbs), len(post_rzbs),[1], [4], ['#FFFFFF'], bccs[10])

        #Build mutation map. This should be the same throughout the EA since SCCs have fixed length
        scc_lengths     = [scc.length for scc in self.sccs]
        self.bcc_to_scc = build_map(scc_lengths, 0)[2]

        #Calculate mutation probs. These will change throughout the EA since the set of valid nucleotide assignments can change
        num_bccs = len ( self.bcc_to_scc.keys() )
        weights  = list()
        for i in range(num_bccs):
            (scc_id, bcc_id_loc) = self.bcc_to_scc[i]
            scc                  = self.sccs[scc_id]
            cur_valid          = scc.valid_nucs[bcc_id_loc]
            if scc.size == 1: # unpaired
                weight = len(cur_valid) - 1  # -1 because we care about the #of valid bccs that a bcc can mutate into 
            else:
                weight = ( len(cur_valid) - 1 ) / float(2) #/2 because each base pair contains 2 bases
            weights.append(weight)
        denom              = sum(weights)
        self.mutation_prob = [w/float(denom) for w in weights]

    def generate_segments(self):
        unordered_segments = list()
        for scc in self.sccs:
            if scc.size == 1: #unpaired
                seq = ''.join(scc.bccs)
                seg = Segment(scc.name[0], scc.length, scc.segs_strand[0], scc.segs_pos[0], scc.segs_colors[0], seq)
                unordered_segments.append(seg)
            elif scc.size == 2:
                seq1 = list()
                seq2 = list()
                for bpp in scc.bccs:
                   seq1.append(bpp[0])
                   seq2.append(bpp[1])
                seq2  = seq2[::-1]
                seq1  = ''.join(seq1)
                seq2  = ''.join(seq2)
                seg1  = Segment(scc.name[0], scc.length, scc.segs_strand[0], scc.segs_pos[0], scc.segs_colors[0], seq1)
                seg2  = Segment(scc.name[1], scc.length, scc.segs_strand[1], scc.segs_pos[1], scc.segs_colors[1], seq2)
                unordered_segments += [seg1, seg2]
                         
        #Sort segments
        self.segments = list()
        self.num_strands = max([seg.strand for seg in unordered_segments]) + 1
        for sd_id in range(self.num_strands):
            strand_segs = [seg for seg in unordered_segments if seg.strand == sd_id]
            sorted_idxs = arg_sort([seg.pos for seg in strand_segs])
            for idx in sorted_idxs:
                self.segments.append( strand_segs[idx] )

        #Extension region seg ids. This is used for novelty distance calculation.
        self.er_seg_ids = list()
        potential_extension_names = ['L0', 'OBS1', 'L1', 'OBS2', 'L2', 'OBS3', 'L3', 'OBS']
        seg_names = [seg.name for seg in self.segments]
        for name in potential_extension_names:
            if name in seg_names:
                idx = inv_list(seg_names, name)
                self.er_seg_ids.append(idx)

    def generate_strands(self):
        self.strands = list()
        for sd_id in range(self.num_strands):
            strand = ''
            strand_segs = [seg for seg in self.segments if seg.strand == sd_id]
            for seg in strand_segs:
                strand += seg.seq
            self.strands.append(strand)

    def update_base_maps(self): #Requires segments to be formed
        lengths                                     = [seg.length for seg in self.segments]
        (ct_base_seg_all, unused, self.base_to_seg) = build_map(lengths, 1)
        self.name_to_bases                          = dict()
        names                                       = [seg.name for seg in self.segments]
        for i, name in enumerate(names):
            self.name_to_bases[name] = ct_base_seg_all[i]

    #This is quite different from the TS version
    def generate_folding_tasks(self, wt_len, stride_size):
        sensor_base_ids = self.name_to_bases['SENSOR'] #These are 1-indexed and right inclusive
        sensor_lims     = (sensor_base_ids[0], sensor_base_ids[-1] + 1) #These are 1-indexed and right exclusive
        sensor_len      = len(sensor_base_ids)
        sd0_len         = len(self.strands[0])
        sd1_len         = len(self.strands[1])

        num_strides  = int((sensor_len - wt_len) / stride_size) + 1

        self.geno = self.strands[0] + "&" +  self.strands[1]
        self.truth_vector = list()
        amp = '&'

        sensor_scc        = self.sccs[4] #!!Hardcoded
        sensor_mismatches = [o + sensor_lims[0] for o in sensor_scc.mismatches]

        logic_tasks = list()
        #no binding
        fc = generate_fold_constraint(sd0_len, [ ]) + amp + generate_fold_constraint(sd1_len, [ ])
        logic_tasks.append((self.geno, fc))
        self.truth_vector.append(0) #We don't want the ribozyme to cleave before the repeats have a chance to bind to the sensor

        #WT binding
        folding_constraints = list()
        lims = [None, None]
        for s in range(num_strides):
            lims[0] = sensor_lims[0] + s*stride_size
            lims[1] = lims[0] + wt_len
            sensor_constraint = [i for i in range(lims[0], lims[1]) if i not in sensor_mismatches] 
            fc = generate_fold_constraint(sd0_len, sensor_constraint ) + amp +  generate_fold_constraint(sd1_len, [ ])
            logic_tasks.append((self.geno, fc))
            self.truth_vector.append(0) #We don't want to cleave the WT

        #Mutant binding
        sensor_constraint = [i for i in range(sensor_lims[0], sensor_lims[1]) if i not in sensor_mismatches] 
        fc = generate_fold_constraint(sd0_len, sensor_constraint ) + amp + generate_fold_constraint(sd1_len, [ ])
        logic_tasks.append((self.geno, fc))
        self.truth_vector.append(1) #We want to cleave the mutant

        self.tasks = [logic_tasks, [], [], []]      
        
    #!!!! TO DO: MOVE THIS SERVER SIDE
    def build_color_string(self):
        self.color_string = ''
        accum = 0 #color string is 0-indexed
        for sd_id in [0,1]:
            strand_segs = [seg for seg in self.segments if seg.strand == sd_id]
            for seg in strand_segs:
                self.color_string += str(accum) + '-' + str(accum + seg.length - 1) + ':' + seg.color + ' ' #?? -1 because color string end is inclusive
                accum += seg.length
            accum += 2 #FORNA requires a gap of 2 between strands in order to display colors correctly

#SCCs, BCCs, valid nucleotide assignments are all explained in "Evolutionary design and analysis of ribozyme‑based logic gates" by N.Kamel
class SCC():
    def __init__(self, name, constraints, length, segs_strand, segs_pos, segs_colors, bccs, mismatches=[]):
        self.name             = name
        self.size             = len(name) #1 for unpaired, 2 for paired
        self.constraints      = constraints
        self.length           = length
        self.segs_strand      = segs_strand #strand on which each segment of the scc is located
        self.segs_pos         = segs_pos    #position within its parent strand of each segment of the scc
        self.segs_colors      = segs_colors
        self.mismatches       = mismatches
        self.bccs             = list()

        self.valid_nucs = list()
        for i in range( self.length ):
            if self.size == 1: #unpaired                
                cur_valid = allowed_lut( self.constraints[i] )
            else:
                bool_mismatch = i in self.mismatches
                cur_valid = allowed_lut( (self.constraints[i][0], self.constraints[i][1], bool_mismatch) )
            self.valid_nucs.append( cur_valid )

        #Executed during initialization only
        if bccs == None:
            for i in range(length):
                self.bccs.append( np.random.choice(self.valid_nucs[i]) )

        #Handle bccs that may become unvalid when an RzBS mutation changes the set of valid assignments
        for i, bcc in enumerate(self.bccs):
            if bcc not in self.valid_nucs[i]:
                self.bccs[i] = np.random.choice(self.valid_nucs[i])
             
class Segment():
    def __init__(self, name, length, strand, pos, color, seq):
        self.name   = name
        self.length = length
        self.strand = strand
        self.pos    = pos
        self.color  = color
        self.seq    = seq
