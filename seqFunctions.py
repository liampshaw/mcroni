from Bio import SeqIO


# Read fasta to dict
def read_fasta(fasta_file):
    '''Reads in fasta file as dict using SeqIO.

    Args:
        fasta_file (str)
            Filename of fasta

    Returns:
        fasta_dict (dict)
            Dictionary of sequences in fasta
    '''
    fasta_parsed = SeqIO.parse(fasta_file, 'fasta') # read in fasta
    fasta_dict = SeqIO.to_dict(fasta_parsed, lambda rec:rec.id) # convert to dict
    return(fasta_dict)


def reverse_complement(seq):
    '''Returns the reverse complement of a DNA sequence.

    Args:
        seq (str)
            String of bases (e.g. ATCCG)

    Returns:
        bases (str)
            Reverse complemented string (e.g. CGGAT)
    '''
    # For reverse_complement function
    alt_map = {'ins':'0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # Replace any insertions.
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases
