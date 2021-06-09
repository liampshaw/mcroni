#!/usr/bin/env python
# vim: set fileencoding=<utf-8> :
# Copyright 2021 Liam Shaw

# Libraries
import sys
import os
from Bio import SeqIO
import subprocess
import re
import argparse
import pandas as pd
import numpy as np
import math

import mcroni.seqFunctions as sf # for conda
# for local usage
#import seqFunctions as sf



def get_options():
    parser = argparse.ArgumentParser(description='Analyse the local genomic context of mcr-1.',
                                     prog='mcroni')
    input_group = parser.add_mutually_exclusive_group(required=True) # mutually exclusive group
    input_group.add_argument('--fasta', help='Fasta file') # either f or l, but not both
    input_group.add_argument('--filelist', help='Alternatively: a list of fasta files')
    parser.add_argument('--output', help='Output file', required=True)
    return parser.parse_args()

def exit_message(message):
    sys.stderr.write(str(message) + "\n")
    sys.exit(1)

# a general function to cut out a section of a genome upstream_bases and downstream_bases away from a gene
def cut_region(fasta_file, upstream_bases=1255, downstream_bases=1867):
    '''Cuts out the region around mcr-1. Default is to cut out the expected size
    of the full composite transposon (e.g. as found in KX528699.1)
    Args:
        fasta_file (str)
            Filename of fasta
        upstream_bases (int)
            Number of bases upstream from *start* of mcr-1. Default=1255
        downstream_bases (int)
            Number of bases downstream from *end* of mcr-1. Default=1867
    Returns:
        region_seq (str)
            Sequence of the requested region
    '''
    if (upstream_bases<0 or downstream_bases<0):
        print('\nWARNING: you specified a negative number of bases. mcroni treats negative bases as 0 (i.e. no flanking region).')
        if upstream_bases<0:
            upstream_bases = 0
        if downstream_bases <0:
            downstream_bases = 0
    print('\nReading in genome from file '+fasta_file+'...')
    contigs = sf.read_fasta(fasta_file)
    print('\nMaking blast database...')
    subprocess.check_call(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl'],\
        stderr=subprocess.DEVNULL,\
        stdout=open(os.devnull, 'w'))
    print('\nSearching for mcr-1...')
    blast_process = subprocess.Popen(['blastn', '-db', fasta_file, \
                            '-query', sf.get_data('mcr1.fa'), \
                            '-outfmt', '6'],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    blast_out, _ = blast_process.communicate() # Read the output from stdout
    blast_output = re.split('\n|\t',blast_out.decode()) # Decode
    print('\nRemoving temporary blast databases...')
    os.remove(fasta_file+'.nin')
    os.remove(fasta_file+'.nhr')
    os.remove(fasta_file+'.nsq')
    # Looking through to cut out upstream region
    if blast_output == ['']:
        print('\nNo blast hit for mcr-1!')
        return
    elif len(blast_output)>13: # should be 12 entries + empty 13th entry for one hit, so more implies >1 hit
        blast_output.pop(-1) # remove empty last entry
        n_hits = int(len(blast_output)/12)
        print('\nWARNING: blast found', n_hits, 'occurrences of mcr-1  in this genome (expected: one hit). This may be biologically interesting but mcroni is not designed for analysing multiple occurrences together. mcroni will analyse ONLY the first hit. Suggest manual inspection and potentially splitting the genome to analyse other occurrences separately.')
    mcr_1_contig = blast_output[1]
    mcr_1_start = int(blast_output[8]) # note that these are base positions so 1-indexed
    mcr_1_end = int(blast_output[9]) # note that these are base positions so 1-indexed
    if int(blast_output[3])==1626:
        if mcr_1_start < mcr_1_end:
            mcr_1_strand = '+' # On positive strand
        else:
            mcr_1_strand = '-' # On negative strand
        print('\nmcr-1 start position is: base', mcr_1_start, 'on the', mcr_1_strand, 'strand', 'of contig', mcr_1_contig)
    else:
        print('\ERROR: sequence does not contain a full-length mcr-1!')
        return
    # Get mcr-1 sequence
    contig_seq = str(contigs[mcr_1_contig].seq)
    if mcr_1_strand=='+':
        mcr_1_seq = contig_seq[mcr_1_start-1:mcr_1_end]
    if mcr_1_strand=='-':
        mcr_1_seq = sf.reverse_complement(contig_seq[mcr_1_end-1:mcr_1_start])
    # POSITIVE STRAND
    if mcr_1_strand == '+':
        upstream_cut_position = (mcr_1_start-1) - upstream_bases # python 0-indexing
        # UPSTREAM
        if upstream_cut_position < 0:
            print('\nWARNING: the mcr-1 contig is not long enough to extract the expected upstream region.')
            print('         --> mcroni will pad the sequence with gaps (-).')
            length_pad = abs(upstream_cut_position)
            mcr_1_upstream = ''.join(['-' for i in range(abs(upstream_cut_position))])+contig_seq[0:(mcr_1_start-1)]
        else:
            mcr_1_upstream = contig_seq[upstream_cut_position:mcr_1_start-1]
        # DOWNSTREAM
        downstream_cut_position = mcr_1_end + downstream_bases
        print('\nCutting out the region from', upstream_cut_position, '-', downstream_cut_position, 'on the positive strand.')
        if downstream_cut_position > len(contig_seq):
            print('\nWARNING: the mcr-1 contig is not long enough to extract the expected downstream region.')
            print('         --> mcroni will pad the sequence with gaps (-).')
            length_pad = downstream_cut_position - len(contig_seq)
            mcr_1_downstream = contig_seq[(mcr_1_end):len(contig_seq)]+''.join(['-' for i in range(abs(length_pad))])
        else:
            mcr_1_downstream = contig_seq[(mcr_1_end):downstream_cut_position]
    # NEGATIVE STRAND
    elif mcr_1_strand == '-':
        # UPSTREAM
        upstream_cut_position = mcr_1_start + upstream_bases
        if upstream_cut_position > len(contig_seq):
            print('\nWARNING: the mcr-1 contig is not long enough to extract the expected upstream region.')
            print('         --> mcroni will pad the sequence with gaps (-).')
            length_pad = upstream_cut_position - len(contig_seq)
            mcr_1_upstream = ''.join(['-' for i in range(length_pad)]) + sf.reverse_complement(contig_seq[mcr_1_start:len(contig_seq)])
        else:
            mcr_1_upstream = sf.reverse_complement(contig_seq[mcr_1_start:upstream_cut_position])
        # DOWNSTREAM
        downstream_cut_position = mcr_1_end -  downstream_bases - 1
        print('\nCutting out the region from', downstream_cut_position, '-', upstream_cut_position, 'on the negative strand.')
        if downstream_cut_position < 0:
            print('\nWARNING: the mcr-1 contig is not long enough to extract the expected downstream region.')
            print('         --> mcroni will pad the sequence with gaps (-).')
            length_pad = abs(downstream_cut_position)
            mcr_1_downstream = sf.reverse_complement(contig_seq[0:mcr_1_end-1])+''.join(['-' for i in range(length_pad)])
        else:
            mcr_1_downstream = sf.reverse_complement(contig_seq[downstream_cut_position:mcr_1_end-1])
    region_seq = mcr_1_upstream+mcr_1_seq+mcr_1_downstream
    return(region_seq)



def cut_upstream_region(fasta_file, threshold=76):
    '''Returns the upstream region of mcr-1 in a genome (assumes just one hit).'''
    print('\nReading in genome from file '+fasta_file+'...')
    contigs = sf.read_fasta(fasta_file)
    print('\nMaking blast database...')
    subprocess.check_call(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl'],\
        stderr=subprocess.DEVNULL,\
        stdout=open(os.devnull, 'w'))
    print('\nSearching for mcr-1...')
    blast_process = subprocess.Popen(['blastn', '-db', fasta_file, \
                            '-query', sf.get_data('mcr1.fa'), \
                            '-outfmt', '6'],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    blast_out, _ = blast_process.communicate() # Read the output from stdout
    blast_output = re.split('\n|\t',blast_out.decode()) # Decode
    print('\nRemoving temporary blast databases...')
    os.remove(fasta_file+'.nin')
    os.remove(fasta_file+'.nhr')
    os.remove(fasta_file+'.nsq')

    # Looking through to cut out upstream region
    if blast_output == ['']:
        print('\nNo blast hit for mcr-1!')
        return
    elif len(blast_output)>13: # should be 12 entries + empty 13th entry for one hit, so more implies >1 hit
        blast_output.pop(-1) # remove empty last entry
        n_hits = int(len(blast_output)/12)
        print('\nWARNING: blast found', n_hits, 'occurrences of mcr-1  in this genome (expected: one hit). This may be biologically interesting but mcroni is not designed for analysing multiple occurrences together. mcroni will analyse ONLY the first hit. Suggest manual inspection and potentially splitting the genome to analyse other occurrences separately.')
    mcr_1_contig = blast_output[1]
    mcr_1_start = int(blast_output[8])
    mcr_1_end = int(blast_output[9])
    if int(blast_output[3])==1626:
        if mcr_1_start < mcr_1_end:
            mcr_1_strand = '+' # On positive strand
        else:
            mcr_1_strand = '-' # On negative strand
        print('\nmcr-1 start position is: base', mcr_1_start, 'on the', mcr_1_strand, 'strand', 'of contig', mcr_1_contig)
    else:
        print('\ERROR: sequence does not contain a full-length mcr-1!')
        return

    # Get mcr-1 variant
    contig_seq = str(contigs[mcr_1_contig].seq)
    if mcr_1_strand=='+':
        mcr_1_seq = contig_seq[mcr_1_start-1:mcr_1_end] # Python indexing
    if mcr_1_strand=='-':
        mcr_1_seq = sf.reverse_complement(contig_seq[mcr_1_end-1:mcr_1_start])
    mcr_1_variant = sf.classify_variant(mcr_1_seq)

    # Use this information to extract the 75bp upstream of mcr-1
    if mcr_1_strand == '+':
        cut_position = mcr_1_start-threshold
        if cut_position < 0:
            print('\nWARNING: the mcr-1 contig is not long enough to extract the expected promoter region.')
            print('nmcroni will pad the sequence with Ns.')
            length_pad = abs(cut_position)
            mcr_1_upstream = ''.join(['N' for i in range(abs(cut_position))])+contig_seq[0:mcr_1_start-1]
        else:
            mcr_1_upstream = contig_seq[cut_position:mcr_1_start-1]


    elif mcr_1_strand == '-':
        cut_position = mcr_1_start+threshold
        if cut_position > len(contig_seq):
            print('\nWARNING: the mcr-1 contig is not long enough to extract the requested upstream region.')
            print('\nmcroni will pad the sequence with Ns.')
            length_pad = cut_position - len(contig_seq) - 1
            mcr_1_upstream = sf.reverse_complement(contig_seq[mcr_1_start:cut_position-1]+''.join(['N' for i in range(length_pad)]))
        else:
            mcr_1_upstream = sf.reverse_complement(contig_seq[mcr_1_start:cut_position-1])
    print('\nThe upstream region of mcr-1 is:')
    print(mcr_1_upstream)
    return([mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, mcr_1_upstream])

def classify_ISApl1_presence(contig, mcr_1_start, mcr_1_strand):
    '''Analyses the upstream and downstream presence of ISApl1 using minimap2.
    Args:
        contig (SeqIO record)
            Record for contig containing mcr-1
        mcr_1_start (int)
            Start position of mcr-1 gene
        mcr_1_strand
            Strand of mcr-1
    Returns:
        ISApl1_dict (dict)
            Dict with keys 'upstream', 'downstream' storing the length, strand of ISApl1
    '''
    # Converting input contig
    if mcr_1_strand == '+':
        contig_str = str(contig.seq)
    if mcr_1_strand == '-': # We construct the sequence so that mcr-1 is on +ve strand
        contig_str = sf.reverse_complement(str(contig.seq))
        mcr_1_start = len(contig_str)-mcr_1_start
    # Write to tmp file
    with open('tmp.fa', 'w') as f:
        f.write('>tmp\n%s' % contig_str)
    # Parameters used in function
    upstream_ISApl1_window = 1264 # Based on 1254 in KX528699
    downstream_ISApl1_window = 3503 # Based on 3493 in KX528699
    positive_map = {True : 'upstream', False : 'downstream'}
    # if ISApl1 strand=='+' then it is right way round, if '-' then opposite way round
    #strand_map = {'+' : 'normal', '-' : 'inverted'}
    orientation_map = {True: 'normal', False: 'inverted'}
    # Gets filled if there are hits
    ISApl1_dict = {'upstream' : [0, 'NA'], 'downstream' : [0, 'NA']}
    # Start positions of ISApl1 will need to be within these limits to count
    upstream_limit = mcr_1_start-upstream_ISApl1_window
    downstream_limit = mcr_1_start+downstream_ISApl1_window
    print('Searching for ISApl1...')
    # Search with blast for ISApl1
    print('\nMaking blast database...')
    subprocess.check_call(['makeblastdb', '-in', 'tmp.fa', '-dbtype', 'nucl'],\
        stderr=subprocess.DEVNULL,\
        stdout=open(os.devnull, 'w'))
    print('\nSearching for ISApl1...')
    blast_process = subprocess.Popen(['blastn', '-db', 'tmp.fa', \
                            '-query', sf.get_data('ISApl1.fa'), \
                            '-outfmt', '6'],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    blast_out, _ = blast_process.communicate() # Read the output from stdout
    blast_output = re.split('\n|\t',blast_out.decode()) # Decode
    print('\nRemoving temporary blast databases...')
    os.remove('tmp.fa'+'.nin')
    os.remove('tmp.fa'+'.nhr')
    os.remove('tmp.fa'+'.nsq')
    print(blast_output)
    # Process blast output
    if blast_output == ['']:
        print('\nNo blast hit for ISApl1.')
        return
    else:
        blast_output.pop(-1) # remove empty trailing entry
        blast_results = pd.DataFrame(np.reshape(blast_output, newshape=(int(np.floor(len(blast_output)/12)), 12)))
        starts = list(pd.to_numeric(blast_results[8])) # will need -1 for python 0-index
        ends = list(pd.to_numeric(blast_results[9])) # will need -1 for python 0-index
        internal_starts = list(pd.to_numeric(blast_results[6])) # internal to ISApl1
        internal_ends = list(pd.to_numeric(blast_results[7])) # internal to ISApl1
        print(list(zip(starts, ends)))
        start_before_end = [starts[i]<ends[i] for i in range(len(starts))]
        orientations = [orientation_map[x] for x in start_before_end] # get the orientation
        print(orientations)
        ISApl1_lengths = [abs(ends[i] - starts[i])+1 for i in range(0, len(ends))] # +1 because e.g. start at 1 finish at 3 means length=3
        print(ISApl1_lengths)
        ISApl1_relative_positions = [positive_map.get(loc, loc) for loc in [x<mcr_1_start for x in ends]]
        ISApl1_limits = [starts[i]>upstream_limit and starts[i]<downstream_limit for i in range(0, len(starts))]
        # Loop through all instances and check if condition is met
        upstream_l = [ISApl1_relative_positions[i]=='upstream' and ISApl1_limits[i] for i in range(0, len(starts))]
        if True in upstream_l:
            upstream_ind = upstream_l.index(True)
            ISApl1_dict['upstream'] = [ISApl1_lengths[upstream_ind],orientations[upstream_ind]]
        downstream_l = [ISApl1_relative_positions[i]=='downstream' and ISApl1_limits[i] for i in range(0, len(starts))]
        if True in downstream_l:
            downstream_ind = downstream_l.index(True)
            ISApl1_dict['downstream'] = [ISApl1_lengths[downstream_ind], orientations[downstream_ind]]
        # return the dict
        print('\nThe summary of ISApl1 presence is:')
        ISApl1_relative_starts = [x-mcr_1_start for x in starts]
        ISApl1_relative_ends = [x-mcr_1_start for x in ends]
        print('Starts:', ISApl1_relative_starts)
        print('Ends:', ISApl1_relative_ends)
        print('Lengths:', ISApl1_lengths)
        presences = list(zip(ISApl1_relative_starts, ISApl1_relative_ends))
        internal_numbers = list(zip(internal_starts, internal_ends))
        print(sorted(presences))
        print([internal_numbers[presences.index(x)] for x in sorted(presences)])
        print(ISApl1_dict)


# Header for output file
output_header = ('\t').join(['FILE', 'ISOLATE', 'MCR1.CONTIG', 'MCR1.START', 'MCR1.STRAND',
                            'MCR1.VARIANT', 'MCR1.UPSTREAM.SEQUENCE',
                            'PLASMIDS.ON.MCR1.CONTIG', 'PLASMIDS.ELSEWHERE',
                            'ISAPL1.UPSTREAM.LENGTH', 'ISAPL1.UPSTREAM.STRAND',
                            'ISAPL1.DOWNSTREAM.LENGTH', 'ISAPL1.DOWNSTREAM.STRAND'])

def main():
    args = get_options()
    print(args)
    fastas = []
    if args.filelist is not None:
        with open(args.filelist, 'r') as f:
            for line in f.readlines():
                fasta_file = line.strip()
                fastas.append(fasta_file)
    else:
        fastas.append(str(args.fasta))


    with open(args.output, 'w') as output_file:
        output_file.write(output_header+'\n')
        for fasta_file in fastas:
            if not os.path.exists(fasta_file):
                print('\nERROR: input fasta does not exist:', fasta_file)
                return
            fasta_name = re.sub('\\..*', '', re.sub('.*\\/', '', fasta_file))
            output = cut_upstream_region(fasta_file)
            if output!=None:
                mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, mcr_1_upstream = output[0], output[1], output[2], output[3], output[4]
                # Write to file
                output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (fasta_file, fasta_name, mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, mcr_1_upstream))
                # Get plasmid replicons too
                plasmids = sf.plasmid_replicons(fasta_file, mcr_1_contig)
                # Write to file
                output_file.write('\t%s\t%s' % (plasmids[0], plasmids[1]))
                # Get ISApl1 status
                ISApl1_status = classify_ISApl1_presence(sf.read_fasta(fasta_file)[mcr_1_contig], mcr_1_start, mcr_1_strand)
                # write to file
                output_file.write('\t%d\t%s\t%d\t%s\n' % (ISApl1_status['upstream'][0], ISApl1_status['upstream'][1], \
                        ISApl1_status['downstream'][0], ISApl1_status['downstream'][1]))

            else:
                output_file.write('%s\t%s\t%s\n' % (fasta_file, fasta_name, '\t'.join(['NA' for i in range(0,11)]))) # 11 empty fields

if __name__ == "__main__":
    main()
