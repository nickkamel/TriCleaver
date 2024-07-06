from .representation import Segment
from ..common.utility import prefix_sum, arg_sort

def generate_segments(sccs, num_strands):
    '''
    Generates segments from segment connected components (SCCs)
    Each SCC specifies either an isolated segment or a pair of reverse complementary segments 

    Parameters:
        sccs (list of SCC objects): The segment connected components
        num_strands (int): Number of strands in the sRz representation (should always be 3)

    Returns:
        segments (list of Segment objects): List of segment objects sorted by strand id first and then by that segment's position within that strand
    '''

    unordered_segments = list()
    for scc in sccs:
        if scc.size == 1: #unpaired
            sequence = ''.join(scc.bccs)
            segment = Segment(scc.name[0], scc.length, scc.strands[0], scc.positions[0], sequence)
            unordered_segments.append(segment)
        elif scc.size == 2:
            sequence_1 = list()
            sequence_2 = list()
            for bcc in scc.bccs:
                sequence_1.append(bcc[0])
                sequence_2.append(bcc[1])
            sequence_2 = sequence_2[::-1]
            sequence_1 = ''.join(sequence_1)
            sequence_2 = ''.join(sequence_2)
            segment_1 = Segment(scc.name[0], scc.length, scc.strands[0], scc.positions[0], sequence_1)
            segment_2 = Segment(scc.name[1], scc.length, scc.strands[1], scc.positions[1], sequence_2)
            unordered_segments += [segment_1, segment_2]
                         
    # Sort segments baesd on their position in their corresponding strand
    segments = list()
    for i in range(num_strands):
        strand_segments = [seg for seg in unordered_segments if seg.strand == i]
        sorted_idxs = arg_sort([seg.pos for seg in strand_segments])
        for idx in sorted_idxs:
            segments.append( strand_segments[idx] )
    return segments


def generate_strands(segments, num_strands):
    '''
    For each strand, concatenates the segment sequences into a single sequence

    Parameters:
        segments (list of Segment objects)
        num_strands (int): Number of strands in the sRz representation (should always be 3 for TC)

    Return:
        strands (list of strings): List containg the sequence of each strand

    '''

    strands = list()
    for sd_id in range(num_strands):
        strand = ''
        strand_segs = [seg for seg in segments if seg.strand == sd_id]
        for seg in strand_segs:
            strand += seg.seq
        strands.append(strand)
    return strands

def build_base_maps(segments):
    '''
    Consider a sequence of sum(lengths) length partitioned into len(lengths) subsequences where the length of subsequence i is lengths[i]
    An element in the sequence has a global position within that sequence as well as a local position within its subsequence
    This function builds maps between these local and global coordinates 

    Parameters: 
        segments (list of Segment objects)

    Returns:
        base_to_seg (dict). Map from integer global nucleotide position to pair (position of the segment within the strand, position of base within that segment) 
        name_to_bases (dict). Map from string segment name to list of nucleotides in that segment in terms of their global posiiton on the strand
    '''

    lengths = [seg.length for seg in segments]
    names = [seg.name for seg in segments]
    num_segs = len(lengths)

    name_to_bases = dict()
    base_to_seg = dict()
    p_sum = [0] + prefix_sum(lengths)
    # This is 1-indexed to interface with the folding output BPPM which is 1-indexed
    p_sum = [p + 1 for p in p_sum]

    for seg_id in range(num_segs):
        seg_ids_global = list(range(p_sum[seg_id], p_sum[seg_id+1]))
        for local_id, global_id in enumerate(seg_ids_global):
            base_to_seg[global_id] = (seg_id, local_id) 
        name_to_bases[names[seg_id]] = seg_ids_global 

    return base_to_seg, name_to_bases
