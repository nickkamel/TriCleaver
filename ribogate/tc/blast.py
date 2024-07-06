import subprocess

def run_blast(sequence, database_name):
    # Prepare the BLAST command
    #cmd = ["blastn", "-query", "-", "-db", database_name]
    #cmd = ["blastn", "-query", "-", "-db", database_name, "-outfmt", "6 qseqid saccver pident qcovs"]  
    cmd = ["blastn", "-query", "-", "-db", database_name, "-outfmt", "6 qseqid saccver pident qcovs", "-task", "blastn-short"] 
    fasta_sequence = f">mySequence\n{sequence}" # Format the sequence in FASTA format

    # Execute the command
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    output, error = process.communicate(input=fasta_sequence)

    # Check for errors
    if process.returncode != 0:
        raise Exception(f"BLAST error: {error}")

    return output

def calculate_specificity(blast_output):
    '''
    Uses the method explained in the supplementary material of Automated design of hammerhead ribozymes and validation by targeting the PABPN1 gene transcript by Nawwaf Kharma, et al
    '''

    specificity_score = 0
    lines = blast_output.split('\n')
    for line in lines:
        if "mySequence" in line:
            cols  = line.split('\t')
            specificity_score += float(cols[2])*float(cols[3])/10000
    return specificity_score