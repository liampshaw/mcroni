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
import logging

# import mcroni.seqFunctions as sf # for conda
# for local usage
import seqFunctions as sf



def get_options():
    parser = argparse.ArgumentParser(description='Analyse the local genomic context of mcr-1.',
                                     prog='mcroni')
    input_group = parser.add_mutually_exclusive_group(required=True) # mutually exclusive group
    output_group = parser.add_mutually_exclusive_group(required=False) # mutually exclusive group
    input_group.add_argument('--fasta', help='Fasta file') # either f or l, but not both
    input_group.add_argument('--filelist', help='Alternatively: a list of fasta files')
    parser.add_argument('--output', help='Output prefix', required=True)
    parser.add_argument('-v', '--verbose', help='verbose output',
                    action='store_true')
    output_group.add_argument('--force', help='Force overwriting of output files.',
                    action='store_true', required=False)
    output_group.add_argument('--append', help='Append to existing output files.',
                                    action='store_true', required=False)
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
        logging.debug('\nWARNING: you specified a negative number of bases. mcroni treats negative bases as 0 (i.e. no flanking region).')
        if upstream_bases<0:
            upstream_bases = 0
        if downstream_bases <0:
            downstream_bases = 0
    logging.debug('\nReading in genome from file '+fasta_file+'...')
    contigs = sf.read_fasta(fasta_file)
    logging.debug('\nMaking blast database...')
    subprocess.check_call(['makeblastdb', '-in', fasta_file, '-dbtype', 'nucl'],\
        stderr=subprocess.DEVNULL,\
        stdout=open(os.devnull, 'w'))
    logging.debug('\nSearching for mcr-1...')
    blast_process = subprocess.Popen(['blastn', '-db', fasta_file, \
                            '-query', sf.get_data('mcr1.fa'), \
                            '-outfmt', '6'],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    blast_out, _ = blast_process.communicate() # Read the output from stdout
    blast_output = re.split('\n|\t',blast_out.decode()) # Decode
    logging.debug('\nRemoving temporary blast databases...')
    os.remove(fasta_file+'.nin')
    os.remove(fasta_file+'.nhr')
    os.remove(fasta_file+'.nsq')
    # Looking through to cut out upstream region
    if blast_output == ['']:
        logging.debug('\nNo blast hit for mcr-1!')
        return
    elif len(blast_output)>13: # should be 12 entries + empty 13th entry for one hit, so more implies >1 hit
        blast_output.pop(-1) # remove empty last entry
        n_hits = int(len(blast_output)/12)
        logging.debug('\nWARNING: blast found', n_hits, 'occurrences of mcr-1  in this genome (expected: one hit). This may be biologically interesting but mcroni is not designed for analysing multiple occurrences together. mcroni will analyse ONLY the first hit. Suggest manual inspection and potentially splitting the genome to analyse other occurrences separately.')
    mcr_1_contig = blast_output[1]
    mcr_1_start = int(blast_output[8]) # note that these are base positions so 1-indexed
    mcr_1_end = int(blast_output[9]) # note that these are base positions so 1-indexed
    if int(blast_output[3])==1626:
        if mcr_1_start < mcr_1_end:
            mcr_1_strand = '+' # On positive strand
        else:
            mcr_1_strand = '-' # On negative strand
        logging.debug('\nmcr-1 start position is: base '+str(mcr_1_start)+' on the '+str(mcr_1_strand)+' strand of contig '+str(mcr_1_contig))
    else:
        logging.debug('\ERROR: sequence does not contain a full-length mcr-1!')
        return
    # Get mcr-1 sequence
    contig_seq = str(contigs[mcr_1_contig].seq)
    if mcr_1_strand=='+':
        mcr_1_seq = contig_seq[mcr_1_start-1:mcr_1_end]
    if mcr_1_strand=='-':
        mcr_1_seq = sf.reverse_complement(contig_seq[mcr_1_end-1:mcr_1_start])
    mcr_1_variant = sf.classify_variant(mcr_1_seq) # classify variant
    # POSITIVE STRAND
    if mcr_1_strand == '+':
        upstream_cut_position = (mcr_1_start-1) - upstream_bases # python 0-indexing
        # UPSTREAM
        if upstream_cut_position < 0:
            logging.debug('\nWARNING: the mcr-1 contig is not long enough to extract the expected upstream region.')
            logging.debug('         --> mcroni will pad the sequence with gaps (-).')
            length_pad = abs(upstream_cut_position)
            mcr_1_upstream = ''.join(['-' for i in range(abs(upstream_cut_position))])+contig_seq[0:(mcr_1_start-1)]
        else:
            mcr_1_upstream = contig_seq[upstream_cut_position:mcr_1_start-1]
        # DOWNSTREAM
        downstream_cut_position = mcr_1_end + downstream_bases
        logging.debug('\nCutting out the region from '+str(upstream_cut_position)+'-'+str(downstream_cut_position)+' on the positive strand.')
        if downstream_cut_position > len(contig_seq):
            logging.debug('\nWARNING: the mcr-1 contig is not long enough to extract the expected downstream region.')
            logging.debug('         --> mcroni will pad the sequence with gaps (-).')
            length_pad = downstream_cut_position - len(contig_seq)
            mcr_1_downstream = contig_seq[(mcr_1_end):len(contig_seq)]+''.join(['-' for i in range(abs(length_pad))])
        else:
            mcr_1_downstream = contig_seq[(mcr_1_end):downstream_cut_position]
    # NEGATIVE STRAND
    elif mcr_1_strand == '-':
        # UPSTREAM
        upstream_cut_position = mcr_1_start + upstream_bases
        if upstream_cut_position > len(contig_seq):
            logging.debug('\nWARNING: the mcr-1 contig is not long enough to extract the expected upstream region.')
            logging.debug('         --> mcroni will pad the sequence with gaps (-).')
            length_pad = upstream_cut_position - len(contig_seq)
            mcr_1_upstream = ''.join(['-' for i in range(length_pad)]) + sf.reverse_complement(contig_seq[mcr_1_start:len(contig_seq)])
        else:
            mcr_1_upstream = sf.reverse_complement(contig_seq[mcr_1_start:upstream_cut_position])
        # DOWNSTREAM
        downstream_cut_position = mcr_1_end -  downstream_bases - 1
        logging.debug('\nCutting out the region from '+str(downstream_cut_position)+'-'+ str(upstream_cut_position)+' on the negative strand.')
        if downstream_cut_position < 0:
            logging.debug('\nWARNING: the mcr-1 contig is not long enough to extract the expected downstream region.')
            logging.debug('         --> mcroni will pad the sequence with gaps (-).')
            length_pad = abs(downstream_cut_position)
            mcr_1_downstream = sf.reverse_complement(contig_seq[0:mcr_1_end-1])+''.join(['-' for i in range(length_pad)])
        else:
            mcr_1_downstream = sf.reverse_complement(contig_seq[downstream_cut_position:mcr_1_end-1])
    region_seq = mcr_1_upstream+mcr_1_seq+mcr_1_downstream
    mcr_1_relative_start = len(mcr_1_upstream)
    # Classifying upstream 'promoter region' variant
    promoter_region = mcr_1_upstream[-75:]
    promoter_variant = sf.classify_variant(promoter_region, sf.get_data('promoter-variants.fa'))
    # Return the results
    return([mcr_1_contig, mcr_1_variant, promoter_variant, mcr_1_start, mcr_1_strand, mcr_1_relative_start, region_seq]) # mcr_1_start is in bases on contig, mcr_1_relative_start is in python for region_seq

def classify_components(region_seq):
    '''Classify the presence of the various components of the mcr-1 composite transposon.
    Args:
        region_seq (str)
            Sequence of the extracted region, such that mcr-1 is on positive strand.
    Returns:
        '''
    return

def classify_ISApl1_presence(region_seq, mcr_1_relative_start):
    '''Analyses the upstream and downstream presence of ISApl1 in the extract region.
    Args:
        region_seq (str)
            Extracted region containing mcr-1.
        mcr_1_relative_start (int)
            Relative start position of mcr-1 gene
    Returns:
        ISApl1_dict (dict)
            Dict with keys 'upstream', 'downstream' storing the length, strand of ISApl1
    '''
    # Write to tmp file
    with open('tmp.fa', 'w') as f:
        f.write('>tmp\n%s' % re.sub('-', '', region_seq)) # sub out the gaps for now
    # Parameters used in function
    upstream_ISApl1_window = 1264 # Based on 1254 in KX528699
    downstream_ISApl1_window = 3503 # Based on 3493 in KX528699
    positive_map = {True : 'upstream', False : 'downstream'}
    # if ISApl1 strand=='+' then it is right way round, if '-' then opposite way round
    #strand_map = {'+' : 'normal', '-' : 'inverted'}
    orientation_map = {True: 'normal', False: 'inverted'}
    # Gets filled if there are hits
    ISApl1_dict = {'upstream': [('NA', 'NA'), ('NA', 'NA')],
                    'downstream': [('NA', 'NA'), ('NA', 'NA')]}
    # Start positions of ISApl1 will need to be within these limits to count
    upstream_limit = mcr_1_relative_start-upstream_ISApl1_window
    downstream_limit = mcr_1_relative_start+downstream_ISApl1_window
    logging.debug('Searching for ISApl1...')
    # Search with blast for ISApl1
    logging.debug('\nMaking blast database...')
    subprocess.check_call(['makeblastdb', '-in', 'tmp.fa', '-dbtype', 'nucl'],\
        stderr=subprocess.DEVNULL,\
        stdout=open(os.devnull, 'w'))
    logging.debug('\nSearching for ISApl1...')
    blast_process = subprocess.Popen(['blastn', '-db', 'tmp.fa', \
                            '-query', sf.get_data('ISApl1.fa'), \
                            '-outfmt', '6'],
                            stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    blast_out, _ = blast_process.communicate() # Read the output from stdout
    blast_output = re.split('\n|\t',blast_out.decode()) # Decode
    logging.debug('\nRemoving temporary blast databases...')
    os.remove('tmp.fa'+'.nin')
    os.remove('tmp.fa'+'.nhr')
    os.remove('tmp.fa'+'.nsq')
    logging.debug(blast_output)
    # Process blast output
    if blast_output == ['']:
        logging.debug('\nNo blast hit for ISApl1.')
        return(ISApl1_dict)
    else:
        blast_output.pop(-1) # remove empty trailing entry
        blast_results = pd.DataFrame(np.reshape(blast_output, newshape=(int(np.floor(len(blast_output)/12)), 12)))
        starts = list(pd.to_numeric(blast_results[8])) # will need -1 for python 0-index
        ends = list(pd.to_numeric(blast_results[9])) # will need -1 for python 0-index
        internal_starts = list(pd.to_numeric(blast_results[6])) # internal to ISApl1
        internal_ends = list(pd.to_numeric(blast_results[7])) # internal to ISApl1
        logging.debug(list(zip(starts, ends)))
        start_before_end = [starts[i]<ends[i] for i in range(len(starts))]
        orientations = [orientation_map[x] for x in start_before_end] # get the orientation
        logging.debug(orientations)
        ISApl1_lengths = [abs(ends[i] - starts[i])+1 for i in range(0, len(ends))] # +1 because e.g. start at 1 finish at 3 means length=3
        logging.debug(ISApl1_lengths)
        ISApl1_relative_positions = [positive_map.get(loc, loc) for loc in [x<mcr_1_relative_start for x in ends]]
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
        logging.debug('\nThe summary of ISApl1 presence is:')
        ISApl1_relative_starts = [x-mcr_1_relative_start for x in starts]
        ISApl1_relative_ends = [x-mcr_1_relative_start for x in ends]
        presences = list(zip(ISApl1_relative_starts, ISApl1_relative_ends))
        logging.debug(presences)
        internal_numbers = list(zip(internal_starts, internal_ends))
        sorted_isapl1_presences = sorted(presences)
        sorted_isapl1_internal = [internal_numbers[presences.index(x)] for x in sorted(presences)]
        logging.debug('\nThe returned dict is:')
        logging.debug(ISApl1_dict)
        ISApl1_summary_list = list(zip(sorted_isapl1_presences, sorted_isapl1_internal))
        if [[x[0]<0 for x in y] for y in ISApl1_summary_list].count([True, False])==1:
            ISApl1_dict['upstream'] = ISApl1_summary_list[0]
            if len(ISApl1_summary_list)>1:
                ISApl1_dict['downstream'] = ISApl1_summary_list[1]
        else:
            ISApl1_dict['downstream'] = ISApl1_summary_list[0]
        return(ISApl1_dict)


# Header for output file
output_header = ('\t').join(['FILE', 'ISOLATE', 'MCR1.CONTIG',
                            'MCR1.START', 'MCR1.STRAND','MCR1.VARIANT', 'PROMOTER.VARIANT',
                            'PLASMIDS.ON.MCR1.CONTIG', 'PLASMIDS.ELSEWHERE',
                            'ISAPL1.UPSTREAM.REL.START', 'ISAPL1.UPSTREAM.REL.END',
                            'ISAPL1.UPSTREAM.INTERNAL.START', 'ISAPL1.UPSTREAM.INTERNAL.END',
                            'ISAPL1.DOWNSTREAM.REL.START', 'ISAPL1.DOWNSTREAM.REL.END',
                            'ISAPL1.DOWNSTREAM.INTERNAL.START', 'ISAPL1.DOWNSTREAM.INTERNAL.END',])

def main():
    args = get_options()
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    fastas = []
    type_of_file_writing = 'w'
    if args.append:
        type_of_file_writing = 'a'

    if args.filelist is not None:
        with open(args.filelist, 'r') as f:
            for line in f.readlines():
                fasta_file = line.strip()
                fastas.append(fasta_file)
    else:
        fastas.append(str(args.fasta))
    output_table_summary = args.output+'_table.tsv'
    output_fasta = args.output+'_sequence.fa'
    if os.path.exists(output_fasta):
        logging.debug('\nWARNING: output fasta already exists')
        if args.force:
            os.remove(output_fasta)
        if args.append:
            logging.debug('\nAppending to existing file.')
        else:
            logging.debug('\nERROR: refusing to overwrite existing output table. If used --force, will append')
            return

    if os.path.exists(output_table_summary):
        logging.debug('\nWARNING: output table already exists.')
        if args.force:
            os.remove(output_table_summary)
        if args.append:
            logging.debug('\nAppending to existing file.')
        else:
            logging.debug('\nERROR: refusing to overwrite existing output table.')
            return
    with open(output_table_summary, type_of_file_writing) as output_file:
        output_file.write(output_header+'\n') # Write the header of table
        for fasta_file in fastas:
            if not os.path.exists(fasta_file):
                logging.debug('\nERROR: input fasta does not exist:', fasta_file)
                return
            fasta_name = re.sub('\\..*', '', re.sub('.*\\/', '', fasta_file))
            output = cut_region(fasta_file)
            if output!=None:
                mcr_1_contig, mcr_1_variant, promoter_variant, mcr_1_start, mcr_1_strand, mcr_1_relative_start, region_seq = output[0], output[1], output[2], output[3], output[4], output[5], output[6]
                # Write to file
                output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (fasta_file, fasta_name, mcr_1_contig, mcr_1_start, mcr_1_strand, mcr_1_variant, promoter_variant))
                # Write sequence too
                if not os.path.exists(output_fasta):
                    with open(output_fasta, 'w') as output_fasta_file:
                        output_fasta_file.write('>%s %s %s\n%s\n' % (fasta_name, mcr_1_contig, mcr_1_variant, region_seq))
                else:
                    with open(output_fasta, 'a') as output_fasta_file:
                        output_fasta_file.write('>%s %s %s\n%s\n' % (fasta_name, mcr_1_contig, mcr_1_variant, region_seq))
                # Get plasmid replicons too
                logging.debug('Finding plasmid replicons using abricate...')
                plasmids = sf.plasmid_replicons(fasta_file, mcr_1_contig)
                # Write to file
                output_file.write('\t%s\t%s' % (plasmids[0], plasmids[1]))
                # Get ISApl1 status
                ISApl1_status = classify_ISApl1_presence(region_seq, mcr_1_relative_start)
                # write to file
                output_file.write('\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (ISApl1_status['upstream'][0][0], ISApl1_status['upstream'][0][1],ISApl1_status['upstream'][1][0], ISApl1_status['upstream'][1][1],
                        ISApl1_status['downstream'][0][0], ISApl1_status['downstream'][0][1], ISApl1_status['downstream'][1][0], ISApl1_status['downstream'][1][1]))
            else:
                output_file.write('%s\t%s\t%s\n' % (fasta_file, fasta_name, '\t'.join(['NA' for i in range(0,14)]))) # 14 empty fields
    if os.path.exists('tmp.fa'): # Clean up tmp file
        os.remove('tmp.fa')

if __name__ == "__main__":
    main()
