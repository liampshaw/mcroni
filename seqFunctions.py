from Bio import SeqIO

# Read fasta to dict
def read_fasta(fasta_file):
    '''Simply uses SeqIO to read in fasta as dict.'''
    return(SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'), lambda rec:rec.id))

# For reverse complementing
alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

def reverse_complement(seq):
    '''Returns the reverse complement of a DNA sequence.'''
    for k,v in alt_map.items():
        seq = seq.replace(k,v)
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    for k,v in alt_map.items():
        bases = bases.replace(v,k)
    return bases
