from folding import generate_fold_constraint, build_map, allowed_lut
import numpy as np
from utility import arg_sort, inv_list

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
        self.fold_constraints    = list()
        self.potential_fitnesses = dict()

        self.folded_strand_ids    = [0, 1] #This is used for building the color string. The ribozyme and the substrate, not the repeats. 

    #In TS each paired SCC explicitly contains the names of both its segments. Here in TC we only give a single name to the SCC and include a procedure here for generating the segment names. 
    def generate_segments(self):
        unordered_segments = list()
        for scc in self.sccs:
            if scc.size == 1: #unpaired
                seq = ''.join(scc.bccs)
                seg = Segment(scc.name, scc.length, scc.segs_strand[0], scc.segs_pos[0], scc.segs_colors[0], seq)
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
                if 'OBS' in scc.name:
                    seg1N = scc.name
                    seg2N = 'Input' + scc.name[-1]
                else: #I'm asssuming the only other paired scc are those of the stems
                    seg1N = scc.name + 'A'
                    seg2N = scc.name + 'B' 
                seg1  = Segment(seg1N, scc.length, scc.segs_strand[0], scc.segs_pos[0], scc.segs_colors[0], seq1)
                seg2  = Segment(seg2N, scc.length, scc.segs_strand[1], scc.segs_pos[1], scc.segs_colors[1], seq2)
                unordered_segments += [seg1, seg2]
                         
        #Sort segments
        self.segments = list()
        self.num_strands = max([seg.strand for seg in unordered_segments]) + 1
        for sd_id in range(self.num_strands):
            strand_segs = [seg for seg in unordered_segments if seg.strand == sd_id]
            sortedIdxs = arg_sort([seg.pos for seg in strand_segs])
            for idx in sortedIdxs:
                self.segments.append( strand_segs[idx] )

        #Extension seg ids
        self.er_seg_ids = list()
        potential_extension_names = ['L0', 'OBS1', 'L1', 'OBS2', 'L2', 'OBS3', 'L3', 'OBS']
        seg_names = [seg.name for seg in self.segments]
        for peName in potential_extension_names:
            if peName in seg_names:
                idx = inv_list(seg_names, peName)
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

    def update_mutation_maps(self, types):
        if 'bccs' in types:
            scc_lengths     = [scc.length for scc in self.sccs]
            self.bcc_to_scc = build_map(scc_lengths, 0)[2]

    def update_mutation_prob(self):
        num_bccs = len ( self.bcc_to_scc.keys() )
        weights  = list()
        for i in range(num_bccs):
            (scc_id, bcc_id_loc) = self.bcc_to_scc[i]
            scc                  = self.sccs[scc_id]
            cur_allowed           = scc.allowed[bcc_id_loc]
            if scc.size == 1: # unpaired
                weight = len(cur_allowed) - 1                # -1 because we care about the #of allowed bccs that a bcc can mutate into 
            else:
                weight = ( len(cur_allowed) - 1 ) / float(2) #/2 because each base pair contains 2 bases
            weights.append(weight)

        denom              = sum(weights)
        self.mutation_prob = [w/float(denom) for w in weights]

    def build_color_string(self):
        self.color_string = ''
        accum = 0 #color string is 0-indexed
        for sd_id in self.folded_strand_ids:
            strand_segs = [seg for seg in self.segments if seg.strand == sd_id]
            for seg in strand_segs:
                self.color_string += str(accum) + '-' + str(accum + seg.length - 1) + ':' + seg.color + ' ' #?? -1 because color string end is inclusive
                accum += seg.length
            accum += 2 #FORNA requires a gap of 2 between strands in order to display colors correctly

    #This is quite different from the TS version
    def generate_folding_tasks(self, wt_len, stride_size):
        obs_base_ids = self.name_to_bases['OBS'] #These are 1-indexed and right inclusive
        obs_lims     = (obs_base_ids[0]   , obs_base_ids[-1] + 1) #These are 1-indexed and right exclusive
        obs_len      = len(obs_base_ids)
        sd0_len      = len(self.strands[0])
        sd1_len      = len(self.strands[1])

        num_strides  = int((obs_len - wt_len) / stride_size) + 1

        self.geno = self.strands[0] + "&" +  self.strands[1]
        self.truth_vector = list()
        amp = '&'

        obs_scc        = self.sccs[4] #!!Hardcoded
        obs_mismatches = [o + obs_lims[0] for o in obs_scc.mismatches]

        logic_tasks = list()
        #no binding
        fc = generate_fold_constraint(sd0_len, [ ]) + amp + generate_fold_constraint(sd1_len, [ ])
        logic_tasks.append((self.geno, fc))
        self.truth_vector.append(0) #We don't want the ribozyme to cleave before the repeats have a change to bind to the OBS

        #WT binding
        fold_constraints = list()
        lims = [None, None]
        for s in range(num_strides):
            lims[0] = obs_lims[0] + s*stride_size
            lims[1] = lims[0] + wt_len
            obs_constraint = [i for i in range(lims[0], lims[1]) if i not in obs_mismatches] 
            fc = generate_fold_constraint(sd0_len, obs_constraint ) + amp +  generate_fold_constraint(sd1_len, [ ])
            logic_tasks.append((self.geno, fc))
            self.truth_vector.append(0) #We don't want to cleave the WT

        #Mutant binding
        obs_constraint = [i for i in range(obs_lims[0], obs_lims[1]) if i not in obs_mismatches] 
        fc = generate_fold_constraint(sd0_len, obs_constraint ) + amp + generate_fold_constraint(sd1_len, [ ])
        logic_tasks.append((self.geno, fc))
        self.truth_vector.append(1) #We want to cleave the mutant

        self.tasks = [logic_tasks, [], [], []]        

class SCC():
    def __init__(self, name, size, spec_constraints, length, segs_strand, segs_pos, segs_colors, seg_lims_arm):
        self.name             = name
        self.size             = size #1 for unpaired, 2 for paired
        self.spec_constraints = spec_constraints #non-wild card constraints
        self.length           = length
        self.segs_strand      = segs_strand #strand on which each segment of the scc is located
        self.segs_pos         = segs_pos    #position within its parent strand of each segment of the scc
        self.segs_colors      = segs_colors
        self.seg_lims_arm     = seg_lims_arm #(start,end) position of segment on arm. None of the segment is not arm. 

        self.constraints      = list()
        self.mismatches       = list()
        self.allowed          = list()
        self.bccs             = list()

    def update_constraints(self, rev_arm_top):
        self.constraints = list()

        #If the scc is not part of the arm, then we store its base constraints in two steps
        #First, we generate a string of wildcard N or (N,N) from that scc's length
        #Then we overwrite it with the specified constraints.
        #Therefore, non-arm scc's with no specified constraints are wildcards.
        #We do things this way to avoid specifying a bunch of wildcards in the template
        #For example, we don't specify any of the inter constraint and we only have to specify the non-wildcard G-C pair for stem 2
        if self.seg_lims_arm == None:
            if self.size == 1: #unpaired
                wild_card_bbp = 'N'
            elif self.size == 2: #paired
                wild_card_bbp = ('N', 'N')
            for i in range(self.length):
                self.constraints.append(wild_card_bbp)

            ##Overwrite Ns with specified constraints (starting from first N)
            for i in range( len(self.spec_constraints) ):
                self.constraints[i] = self.spec_constraints[i] 
        else:
            (a, b) = (self.seg_lims_arm[0], self.seg_lims_arm[1])
            if self.size == 2:
                self.constraints    = [('N', base) for base in rev_arm_top[a:b] ]
                #if (a, b) == (8, 13): #!!! Added May 2020. Hard coded closing stem 3 A constraint. stem 3
                if self.name == 'stem3':
                    self.constraints[0] = ('A', self.constraints[0][1]) #This should always be ('A', 'U') 
            else:
                self.constraints = [base for base in rev_arm_top[a:b] ]

    def update_allowed(self): #(during initialization, size change, mismatch change)
        self.allowed = list()
        for i in range( self.length ):
            if self.size == 1: #unpaired                
                cur_allowed = allowed_lut( self.constraints[i] )
            else:
                bool_mismatch = i in self.mismatches
                cur_allowed = allowed_lut( (self.constraints[i][0], self.constraints[i][1], bool_mismatch) )
            self.allowed.append( cur_allowed )


        #Handle length constraint violations
        num_bccs = len(self.bccs)
        new_len  = len(self.allowed)
        if num_bccs < new_len:
            for i in range(num_bccs, new_len):
                self.bccs.append( np.random.choice(self.allowed[i]) )
        if num_bccs > new_len:
            self.bccs = self.bccs[0:new_len]

        #Handle bbp constraint violations
        for i, bcc in enumerate(self.bccs):
            if bcc not in self.allowed[i]:
                self.bccs[i] = np.random.choice(self.allowed[i])
               

class Segment():
    def __init__(self, name, length, strand, pos, color, seq):
        self.name   = name
        self.length = length
        self.strand = strand
        self.pos    = pos
        self.color  = color
        self.seq    = seq
