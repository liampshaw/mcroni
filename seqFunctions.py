#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

# Useful functions for DNA sequences

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
    '''Returns the reverse complement of a DNA sequence. Assumes no insertions.

    Args:
        seq (str)
            String of bases (e.g. ATCCG)

    Returns:
        bases (str)
            Reverse complemented string (e.g. CGGAT)
            (N.B. if seq contains non-standard characters, returns None)
    '''
    # Mapping of bases to complement (note that N->N)
    complement_map = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    bases = list(seq.upper()) # Make sure upper-case
    # Check if non-standard bases are present, and error if so
    if any([x not in list(complement_map.keys()) for x in bases]):
        non_standard_bases = [x for x in bases if x not in list(complement_map.keys())]
        print('Refusing to reverse complement!')
        print('Your sequence contains non-standard characters:', ''.join(set(non_standard_bases)))
        return
    else:
        bases = reversed([complement_map.get(base,base) for base in bases])
        bases = ''.join(bases)
        return(bases)
