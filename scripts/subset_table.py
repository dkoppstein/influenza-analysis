'''
This file takes the "full table" of mapped influenza sequences and metainformation, and subsets it based on various criteria. 
'''

import re
import pandas as pd
import argparse
import sys
import matplotlib.pyplot as plt
from mpltools import style
from Bio.SeqRecord import SeqRecord
from pssm_plot import explode_series, nucleotide_frequencies
from cutadapt.align import globalalign_locate
from scipy.stats import rankdata
import numpy as np
from graph_of_bedtools_intersect import calculate_num_times_map

def seqs_align(seq1, seq2, error_rate=0.1):
    '''Do sequences 1 and 2 align given the following criteria:

    1. Error rate = 0.1 per 10 nucleotides (int(floor(0.1 * len(seq))))
    2. The alignments perfectly match up
    '''
    len_seq1 = len(seq1)
    if len_seq1 != len(seq2):
        return False
    # use C bindings for superfast alignment!
    aln = globalalign_locate(seq1, seq2, error_rate) 
    # if start1 = start2 and stop1 = stop2 and
    return aln[0] == aln[2] and aln[1] == aln[3]

def restrict_table_by_length(table, l):
    return table[table['ORIGINAL_SEQUENCE_LENGTH'] >= l]

def limit_num_maps(table, n):
    'Returns a table with only the reads that map <= n times to TSSes'
    return table.groupby('READ_NAME').filter(lambda x: len(x) <= n)

def genes_are_same_snrna(first, second):
    'Deals with gene names that are something like RNU7-12P, etc.'
    regex = 'U\d[\d]?'
    first, second = tuple(re.search(regex, n) for n in (first, second))
    return first is not None and second is not None and \
      first.group() == second.group()

def drop_paralogs(df, error_rate=0.2):
    '''Given a DataFrame (preferably grouped by read_name,
    drops paralogous rows, as determined by

    a) sequences threeprime of cleavage aligning (as determined by seqs_align), and
    b) the gene types are the same.
    c) gene_name is the same
    
    O(n^2) operation. NB: does this in place!'''
    i, j = 0, 1
    while i < len(df):
        while j < len(df):
            i_lab, j_lab = df.index[i], df.index[j]
            i_seq, j_seq = (df['THREEPRIME_OF_CLEAVAGE'][i_lab],
                            df['THREEPRIME_OF_CLEAVAGE'][j_lab])
            i_gn, j_gn = (df['GENE_NAME'][i_lab],
                          df['GENE_NAME'][j_lab])
            if (i_gn == j_gn) or genes_are_same_snrna(i_gn, j_gn) or \
              (df['GENE_TYPE'][i_lab] == df['GENE_TYPE'][j_lab] and \
                i_seq is not np.nan and j_seq is not np.nan and \
                seqs_align(i_seq, j_seq, error_rate=error_rate)):
                df.drop(j_lab, axis=0, inplace=True)
            else:
                j += 1
        i += 1
        j = i + 1
    df.index = range(len(df))
    return df

def drop_all_paralogs(df):
    '''This applies the previous function to each mapping grouped by read_name,
    so as not to lose the associated metainformation of each read.'''
    return df.groupby('READ_NAME').apply(lambda x: drop_paralogs(x))

def process_file(args):
    df = pd.read_table(args.infile, index_col=None)
    if args.drop_paralogs:
        df = drop_all_paralogs(df)
    if args.min_length:
        df = df[df['ORIGINAL_SEQUENCE_LENGTH'] >= args.min_length]
    if args.restrict_to_zero:
        df = df[df['OFFSET_FROM_START'] == 0]
    if args.restrict_around:
        df = df[(df['OFFSET_FROM_START'] >= args.restrict_around[0]) & \
                (df['OFFSET_FROM_START'] <= args.restrict_around[1])]
    if args.remove_sequence:
        df = df[df['ORIGINAL_SEQUENCE'] != args.remove_sequence]
    if args.endswith:
        criterion = df['ORIGINAL_SEQUENCE'].map(
            lambda x: x.endswith(args.endswith))
        df = df[criterion]
    if args.threeprime_startswith:
        df = df[df['THREEPRIME_OF_CLEAVAGE'].map(
            lambda x: x.startswith(args.threeprime_startswith))]
    if args.three_prime_trimmed is not None:
        if args.three_prime_trimmed == '':
            df = df[df['3PTRIMMED'].isnull()]
        else:
            df = df[df['3PTRIMMED'] == args.three_prime_trimmed]
    if args.normalize_by_num_maps:
        df = calculate_num_times_map(df)
        df['WEIGHT'] = df['COUNT'] / df['NUM_TIMES_MAP']
        df['WEIGHT'] = df['WEIGHT'].fillna(0)
    if args.use_rank:
        if 'WEIGHT' in df.columns:
            df['WEIGHT'] = rankdata(df['WEIGHT'])
        else:
            df['WEIGHT'] = rankdata(df['COUNT'])
    if not ('WEIGHT' in df.columns):
        df['WEIGHT'] = df['COUNT']
    # combine upstream and downstream into one sequence
    if args.combine_sequences:
        left = df['ORIGINAL_SEQUENCE'].map(lambda x: x[args.combine_sequences[0]:])
        right = df['THREEPRIME_OF_CLEAVAGE'].map(lambda x: x[:args.combine_sequences[1]])
        df['FINAL_SEQUENCE']  = left + right
    # get the selected columns
    if 'all' not in args.columns:
        df = df[args.columns]

    # get rid of empty and null sequences
    if args.threeprime_of_cleavage:
        nucs = df[['THREEPRIME_OF_CLEAVAGE', 'COUNT']]
        nucs = nucs[nucs['THREEPRIME_OF_CLEAVAGE'] != '']
        nucs = nucs[not nucs['THREEPRIME_OF_CLEAVAGE'].isnull()]
        export = explode_series(nucs['THREEPRIME_OF_CLEAVAGE'])
        export = nucleotide_frequencies(export, count_series=nucs['COUNT'],
                                        normalize=True, ignore_ns=args.ignore_ns)
        export.to_csv(args.outfile, sep='\t')
    else:
        df.to_csv(args.outfile, sep='\t', index=False)
    # if no_modulate_by_count, count is None so get_nucleotide_frequencies
    # does not modulate by count_series

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-l', '--min-length', type=int, nargs='?',
                        required=False)
    parser.add_argument('--threeprime_startswith', default=None, type=str)
    parser.add_argument('--restrict-around', nargs=2, default=None, type=int)
    parser.add_argument('--restrict-to-zero', action='store_true',
                        help=('If specified, will only count those that map to a'
                        'zero position on the TSS.'))
    parser.add_argument('--remove-sequence', default=None, type=str)
    parser.add_argument('--endswith', type=str, nargs='?')
    parser.add_argument('--use-rank', action='store_true')
    parser.add_argument('--drop-paralogs', action='store_true')
    parser.add_argument('--three-prime-trimmed', type=str, nargs='?',
                        required=False)
    parser.add_argument('--threeprime-of-cleavage', action='store_true')
    parser.add_argument('--normalize-by-num-maps', action='store_true')
    parser.add_argument('--as-pssm', action='store_true')
    parser.add_argument('--ignore-ns', action='store_true')
    parser.add_argument('--combine-sequences', default=None, nargs=2, type=int)
    parser.add_argument('--columns', type=str, nargs='+', default=['all'],
                        help=('If want all columns, specify \'all\''))
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
