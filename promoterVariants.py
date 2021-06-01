#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

# Libraries
import sys
import os
from Bio import SeqIO
import subprocess
import re
import pandas as pd
import numpy as np

import seqFunctions as sf

def plasmid_replicons(fasta_file, contig_name, database='plasmidfinder'):
    '''Identifies and returns a list of plasmid replicons present on a contig (using abricate) and also within the wider genome.

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

def classify_mcr_1_variant(mcr_1_seq):
    '''Classifies an mcr-1 variant (mcr-1.1 to mcr-1.13 included).

    Args:
        mcr_1_seq (str)
            Sequence of mcr-1 gene (must be full-length).

    Returns:
        variant_name (str)
            Name of variant. 'Other' if sequence non-identical to any known variant.
    '''
    mcr_1_variants = SeqIO.to_dict(SeqIO.parse('mcr-1-variants.fasta', 'fasta'),
                                lambda rec : rec.id)
    for variant in mcr_1_variants:
        if mcr_1_seq==str(mcr_1_variants[variant].seq):
            variant_name = re.sub('\\|.*', '', variant)
        else:
            variant_name = 'Other' # If sequence isn't a named variant
    return(variant_name)

def classify_ISApl1_presence(contig, mcr_1_start, mcr_1_strand):
    '''Analyses the upstream and downstream presence of ISApl1 using minimap2.
    Returns a dict with the length of ISApl1 and the strand.
    Involves a collapsing of information - doesn't return start and end positions.'''
    # Search with minimap2 for ISApl1
    ISApl1_dict = {'upstream' : [0, 'NA'], 'downstream' : [0, 'NA']}
    if mcr_1_strand == '+':
        contig_str = str(contig.seq)
    if mcr_1_strand == '-':
        contig_str = reverse_complement(str(contig.seq))
        mcr_1_start = len(contig_str)-mcr_1_start
    # Need two tests for ISApl1 - upstream of mcr-1? if yes, then how much? then, downstream of mcr-1, and if yes, then how much?
    upstream_ISApl1_window = 1260 # 1254 in KX528699
    downstream_ISApl1_window = 3500 # 3493 in KX528699
    print('Searching for ISApl1...')
    f = open('tmp.fa', 'w')
    f.write('>tmp\n%s' % contig_str)
    f.close()
    minimap_process = subprocess.Popen(['minimap2', 'ISApl1.fasta', \
                            'tmp.fa'],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    minimap_out, _ = minimap_process.communicate() # Read the output from stdout
    minimap_output = re.split('\t|\n', minimap_out.decode()) # decode
    minimap_output.remove('')
    os.remove('tmp.fa') # remove file
    if minimap_output!=[]:
        # Minimap output has 18 entries. Split based on this
        minimap_results = pd.DataFrame(np.reshape(minimap_output, newshape=(int(np.floor(len(minimap_output)/18)), 18)))
        starts = list(pd.to_numeric(minimap_results[2]))
        ends = list(pd.to_numeric(minimap_results[3]))
        strands = list(minimap_results[4])
        # name seq_length start end strand
        # get lengths of ISApl1 hits
        ISApl1_lengths = [abs(ends[i] - starts[i]) for i in range(0, len(ends))] # check +1 etc for precise length of hits. abs takes care of strand
        # get relative position to mcr-1 - could need to check circularity...although this will rarely be a problem it could conceivably happen
        # Map greater than mcr-1 start position context to
        # Use the relative end of ISApl1 compared to mcr-1 to find out if a hit is upstream or downstream
        positive_map = {True : 'upstream', False : 'downstream'}
        ISApl1_relative_positions = [positive_map.get(loc, loc) for loc in [x<mcr_1_start for x in ends]]

        # Start positions need to be within these limits to count
        upstream_limit = mcr_1_start-upstream_ISApl1_window
        downstream_limit = mcr_1_start+downstream_ISApl1_window
        ISApl1_limits = [starts[i]>upstream_limit and starts[i]<downstream_limit for i in range(0, len(starts))]
        # Loop through all instances and check if condition is met
        upstream_l = [ISApl1_relative_positions[i]=='upstream' and ISApl1_limits[i] for i in range(0, len(starts))]
        # if ISApl1 strand=='+' then it is right way round, if '-' then opposite way round
        strand_map = {'+' : 'normal', '-' : 'inverted'}
        if True in upstream_l:
            upstream_ind = upstream_l.index(True)
            ISApl1_dict['upstream'] = [ISApl1_lengths[upstream_ind],strand_map[strands[upstream_ind]]]
        downstream_l = [ISApl1_relative_positions[i]=='downstream' and ISApl1_limits[i] for i in range(0, len(starts))]
        if True in downstream_l:
            downstream_ind = downstream_l.index(True)
            ISApl1_dict['downstream'] = [ISApl1_lengths[downstream_ind], strand_map[strands[downstream_ind]]]
    # return the dict
    return(ISApl1_dict)



def cut_upstream_region(fasta_file, threshold=76):
    '''Returns the upstream region of mcr-1 in a genome (assumes just one hit).'''
    print('Reading in genome from file '+fasta_file+'...')
    contigs = sf.read_fasta(fasta_file)
    print('Making blast database...')
    subprocess.check_call(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl'],\
        stderr=subprocess.DEVNULL,\
        stdout=open(os.devnull, 'w'))
    print('Searching for mcr-1...')
    blast_process = subprocess.Popen(['blastn', '-db', fasta_file, \
                            '-query', 'mcr-1.fasta', \
                            '-outfmt', '6'],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    blast_out, _ = blast_process.communicate() # Read the output from stdout
    blast_output = blast_out.decode().split('\t') # Decode
    print(blast_output)
    print('Removing temporary blast databases...')
    os.remove(fasta_file+'.nin')
    os.remove(fasta_file+'.nhr')
    os.remove(fasta_file+'.nsq')

    # Looking through to cut out upstream region
    if blast_output == ['']:
        print('No blast hit for mcr-1!')
        return
    mcr_1_contig = blast_output[1]
    mcr_1_start = int(blast_output[8])
    mcr_1_end = int(blast_output[9])
    if int(blast_output[3])==1626:
        if mcr_1_start < mcr_1_end:
            mcr_1_strand = '+' # On positive strand
        else:
            mcr_1_strand = '-' # On negative strand
        print('mcr-1 start position is:', mcr_1_start, 'on the', mcr_1_strand, 'strand', 'of contig', mcr_1_contig)
    else:
        print('Does not contain a full-length mcr-1!')
        return

    # Get mcr-1 variant
    contig_seq = str(contigs[mcr_1_contig].seq)
    if mcr_1_strand=='+':
        mcr_1_seq = contig_seq[mcr_1_start-1:mcr_1_end]
    if mcr_1_strand=='-':
        mcr_1_seq = reverse_complement(contig_seq[mcr_1_end-1:mcr_1_start])
    mcr_1_variant = classify_mcr_1_variant(mcr_1_seq)

    # Use this information to extract the 75bp upstream of mcr-1
    if mcr_1_strand == '+':
        cut_position = mcr_1_start-threshold
        if cut_position < 0:
            print('Contig is not long enough...')
            return([mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, ''])
        mcr_1_upstream = contig_seq[cut_position:mcr_1_start-1]


    elif mcr_1_strand == '-':
        cut_position = mcr_1_start+threshold
        if cut_position > len(contig_seq):
            print('Contig is not long enough...')
            return([mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, ''])
        mcr_1_upstream = reverse_complement(contig_seq[mcr_1_start:cut_position-1])
    print(mcr_1_upstream)
    return([mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, mcr_1_upstream])

if __name__ == "__main__":
    # First argument is input filename (fasta)
    # Second argument is output filename ()
    with args[1]

with open(sys.argv[1], 'r') as f:
    with open(sys.argv[2], 'w') as output_file:
        output_header = 'file\tsample\tcontig\tmcr1.start\tmcr1.strand\tmcr1.variant\tmcr1.upstream.seq\tplasmids.contig\tplasmids.elsewhere\tisapl1.upstream.length\tisapl1.upstream.orientation\tisapl1.downstream.length\tisapl1.downstream.length\n'
        output_file.write(output_header)
        for line in f.readlines():
            fasta_file = line.strip()
            fasta_name = re.sub('\\..*', '', re.sub('.*\\/', '', fasta_file))
            output = cut_upstream_region(fasta_file)
            if output!=None:
                mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, mcr_1_upstream = output[0], output[1], output[2], output[3], output[4]
                # Write to file
                output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (fasta_file, fasta_name, mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, mcr_1_upstream))
                # Get plasmid replicons too
                plasmids = plasmid_replicons(fasta_file, mcr_1_contig)
                # Write to file
                output_file.write('\t%s\t%s' % (plasmids[0], plasmids[1]))
                # Get ISApl1 status
                ISApl1_status = classify_ISApl1_presence(sf.read_fasta(fasta_file)[mcr_1_contig], mcr_1_start, mcr_1_strand)
                # write to file
                output_file.write('\t%d\t%s\t%d\t%s\n' % (ISApl1_status['upstream'][0], ISApl1_status['upstream'][1], \
                        ISApl1_status['downstream'][0], ISApl1_status['downstream'][1]))

            else:
                output_file.write('%s\t%s\t%s\n' % (fasta_file, fasta_name, '\t'.join(['NA' for i in range(0,11)]))) # 11 empty fields
