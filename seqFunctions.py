#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

# Useful functions for DNA sequences

from Bio import SeqIO
import os
import subprocess
import pandas as pd

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
    path_to_file = os.path.abspath(fasta_file) # get absolute path
    fasta_parsed = SeqIO.parse(path_to_file, 'fasta') # read in fasta
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

def classify_variant(seq, variants_db='data/mcr1-variants.fa'):
    '''Classifies a variant sequence from a db fasta of known variants.

    Args:
        seq (str)
            Sequence of the gene/region (must be full-length).
        variants_db (str)
            Fasta filename of known variant database.
                mcr1-variants.fa (default)
                    mcr-1 gene. mcr-1.1 to mcr-1.13 included.
                promoter-region-variants.fa
                    75 bp region upstream of mcr-1 gene.
                    Consensus sequence + 8 variants as classified by
                    Lois Ogunlana, University of Oxford.

    Returns:
        variant_name (str)
            Name of variant. 'Other' if sequence non-identical to any known variant.
                - to add: return closest variant(s)
    '''
    variants = SeqIO.to_dict(SeqIO.parse(variants_db, 'fasta'),
                                lambda rec : rec.id)
    for variant in variants:
        if seq==str(variants[variant].seq):
            variant_name = re.sub('\\|.*', '', variant)
        else:
            variant_name = 'Other' # If sequence isn't a named variant
    return(variant_name)

def plasmid_replicons(fasta_file, contig_name, database='plasmidfinder'):
    '''Identifies and returns a list of plasmid replicons present on a contig and also within the whole fasta.

    Args:
        fasta_file (str)
            Filename of fasta
        contig_name (str)
            ID of contig within fasta
        database (str)
            Abricate database to use

    Returns:
        plasmid_output (list)
            List of two comma-separated strings:
                - plasmid replicons found on contig_name
                - plasmid replicons found elsewhere in genome

    N.B. Uses default thresholds of abricate (--minid 80, --mincov 80).
    '''
    print('Finding plasmid replicons using abricate...')
    abricate_command = ['abricate', '--db', database, fasta_file]
    abricate_process = subprocess.Popen(abricate_command, stdout = subprocess.PIPE, stderr = subprocess.PIPE) # Runs abricate
    abricate_out, _ = abricate_process.communicate() # Read the output from stdout
    abricate_out = abricate_out.decode() # Decode
    with open('tmp.out', 'w') as out: # Write to file
        out.write(abricate_out)
    # Read back in
    data = pd.read_csv('tmp.out', sep='\t', header = 0)
    # Get all plasmid replicons elsewhere in genome
    total_replicons = list(data[~data['SEQUENCE'].str.contains(contig_name, regex=False)]['GENE'])
    # Get just those on the contig
    contig_replicons = list(data[data['SEQUENCE'].str.contains(contig_name, regex=False)]['GENE'])
    os.remove('tmp.out') # Remove abricate file
    # Output list of two comma-separated strings
    plasmid_output = [','.join(contig_replicons), ','.join(total_replicons)]
    return(plasmid_output)
