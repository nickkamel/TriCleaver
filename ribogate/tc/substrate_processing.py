def extract_every_rzbs(downstream_seq, s3_len, rzbs_len, rzbs_extension_len):
    rzbs_positions = []
    rzbs_cut_site_motifs = []
    for i in range(rzbs_extension_len, len(downstream_seq) - rzbs_len  - rzbs_extension_len): # Makes sure the rzbs extension always fits in downstream seq
        cut_site_motif = downstream_seq[i + s3_len - 2:i + s3_len + 1]
        if cut_site_motif[1] == 'U' and cut_site_motif[2] != 'G': # NUH
            rzbs_positions.append(i)
            rzbs_cut_site_motifs.append(cut_site_motif)
    return rzbs_positions, rzbs_cut_site_motifs

def count_tandem_repeats(sequence, repeat_unit_len=3):
    '''
    Identifies the most common tandem repeat in a sequence and returns its count. Considers every possible frame shift.
    
    Parameters:
        sequence (string): The sequence in which the repeats are located
        repeat_unit_len: The length of the repeat units we are counted. For example, use 3 to trinucleotide repeats.
        
    Returns:
        most_common_unit (string): Sequence of the most common repeat unit
        most_common_unit_count (int): Number of instances of the most common unit.
    '''

    # When counting the number of repeats for a repeat unit of length N, we need to consider N different frames 
    frame_offsets = list(range(repeat_unit_len)) # i.e. for triplets we have [0,1,2]
    most_common_unit = None
    most_common_unit_count = 0
    for offset in frame_offsets:
        frame_units = [sequence[i:i+repeat_unit_len] for i in range(offset, len(sequence), repeat_unit_len)]
        frame_units = [t for t in frame_units if len(t) == repeat_unit_len] # Handle case where the len(sequence) - offset is not a mulitiple of repeat_unit_len
        unit_counts = {}
        for unit in frame_units:
            if unit not in unit_counts:
                unit_counts[unit] = 0
            unit_counts[unit] += 1
        # Identify the most common unit within the current frame
        most_common_frame_unit = max(unit_counts, key=unit_counts.get)
        
        # Determine whether the most common unit of these frame is more common than the units of the previous frames
        if unit_counts[most_common_frame_unit] > most_common_unit_count:
            most_common_unit = most_common_frame_unit
            most_common_unit_count = unit_counts[most_common_frame_unit]
    
    return (most_common_unit, most_common_unit_count)
