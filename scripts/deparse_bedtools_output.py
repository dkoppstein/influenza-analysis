from __future__ import absolute_import

import sys
import pandas as pd
import flulib
from Bio import SeqIO
import argparse
import numpy as np
import time
if pd.__version__ < '0.13.1':
    raise ImportError('Requires pandas 0.13.1.')
from cutadapt.align import globalalign_locate
from math import floor

'''
This file takes bedtools output and creates a single table with all
possible information about a set of reads.
'''

BEDTOOLS_INTERSECT_COLUMNS = ['READ_CHROM', 'READ_START', 'READ_END',
                              'READ_NAME', 'READ_CODE', 'READ_STRAND',
                              'TSS_CHROM', 'TSS_START', 'TSS_END', 'TSS_NAME',
                              'TSS_CODE', 'TSS_STRAND', 'TSS_SOURCE',
                              'TSS_ATTRIBUTE', 'UNKNOWN_DOT', 'TSS_META']

COLUMNS_TO_ADD = ['ORIGINAL_SEQUENCE', 'ORIGINAL_SEQUENCE_LENGTH',
                  'OFFSET_FROM_START', '5PTRIMMED', '3PTRIMMED',
                  'GENOMIC_SEQUENCE', 'COUNT', 'SEQUENCE_COUNT',
                'ZEROS', 'ZEROS_TO_TOTAL', 'NUM_READS_PER_SEQUENCE',
                'THREEPRIME_OF_CLEAVAGE']

GTF_KEYS = ['GENE_STATUS', 'EXON_NUMBER', 'LEVEL', 'TRANSCRIPT_TYPE', 'TAG',
            'GENE_ID', 'TRANSCRIPT_ID', 'HAVANA_TRANSCRIPT', 'HAVANA_GENE',
            'TRANSCRIPT_NAME', 'GENE_TYPE', 'TRANSCRIPT_STATUS', 'GENE_NAME']

def num_reads_per_sequence(table):
    '''
    Returns a Series of the number of reads each sequence is represented by.
    sequence -> num_reads
    '''
    groupby_seq = table.groupby('ORIGINAL_SEQUENCE')
    d = dict()
    for original_sequence, group in groupby_seq:
        size = group.groupby('READ_NAME').size()
        d[original_sequence] = size
    return pd.Series(d)
    
def get_sequence_abundances(table):
    '''
    Returns the abundance of each sequence, in counts, as a pandas Series
    sequence -> count
    '''
    no_dups = table.drop_duplicates(cols='READ_NAME', take_last=True)
    return no_dups['COUNT'].astype(int).groupby(no_dups['ORIGINAL_SEQUENCE']).sum()

def dict_of_seqs(f):
    '''Gives a dictionary of the sequences in the given fasta file,
    where the key is a tuple of the chrom, start, stop, strand'''
    d = dict()
    for seq_record in SeqIO.parse(f, format='fasta'):
        tup = flulib.parse_bowtie_interval(seq_record.id)
        d[tup] = str(seq_record.seq).upper()
    return d

def get_sequence(row, d, prefix='TSS'):
    '''Gets a sequence from dictionary d using the information in
    row "row" prefixed by "prefix".'''
    tuple_fields = (prefix + s for s in ('_CHROM', '_START',
                                         '_END', '_STRAND'))
    tup = tuple(str(row[x]) for x in tuple_fields)
    return d[tup]

def get_threeprime_sequence(row, genome_dict, num_nucs=50):
    chrom = row['TSS_CHROM']
    if chrom == '.':
        return None
    if row['READ_STRAND'] == '+':
        end = int(row['READ_END'])
        return str(genome_dict[chrom][end:end+num_nucs].seq).upper()
    else:
        end = int(row['READ_START'])
        return str(genome_dict[chrom][end-num_nucs:end].reverse_complement().seq).upper()

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
            if (i_gn == j_gn) or (df['GENE_TYPE'][i_lab] == df['GENE_TYPE'][j_lab] and \
                i_seq is not np.nan and j_seq is not np.nan and \
                seqs_align(i_seq, j_seq, error_rate=error_rate)):
                df.drop(j_lab, axis=0, inplace=True)
            else:
                j += 1
        i += 1
        j = i + 1
    df.index = range(len(df))
    return df

def keep_only_with_g_after(df):
    '''Discards reads that don't have a G after them when mapped.'''

def drop_all_paralogs(df):
    '''This applies the previous function to each mapping grouped by read_name,
    so as not to lose the associated metainformation of each read.'''
    return df.groupby('READ_NAME').apply(lambda x: drop_paralogs(x))

def remove_duplicate_mappings(table):
    '''
    NB: DEPRECATED IN FAVOR OF DROP_ALL_PARALOGS! 

    For each read, removes mapping events that have duplicate "three prime of
    cleavage" sequences. In other words, if a read maps to multiple U3
    paralogs, then it will get rid of these.
    '''
    grouped = table.groupby('READ_NAME')
    df = grouped.apply(lambda x: x.drop_duplicates(cols='THREEPRIME_OF_CLEAVAGE'))
    df.index = range(len(df)) # reset the index
    return df
    
def process_files(args):
    table = pd.read_table(args.infile, index_col=None)
    # Make dictionary of TSS sequences
    genome_dict = SeqIO.to_dict(SeqIO.parse(args.fasta_file[0], 'fasta'))
    original_reads_dict = SeqIO.index(args.original_reads_file[0], 'fastq-illumina')

    # The number of times a read maps in total
    num_times_map = table.groupby('READ_NAME').size()
    
    # index by read name
    table = table.set_index('READ_NAME')
    
    # original sequence (keyed by read name)
    def get_seq_from_original_dict(row):
        return str(original_reads_dict[row.name].seq)
    
    table['ORIGINAL_SEQUENCE'] = table.apply(get_seq_from_original_dict, axis=1)
    table['ORIGINAL_SEQUENCE_LENGTH'] = table['ORIGINAL_SEQUENCE'].map(len)

    # get dist to mid
    def offset_from_start(row):
        dist_to_mid = (row['TSS_END'] - row['TSS_START']) / 2
        if row['READ_STRAND'] == '+':
            tss_pos = row['TSS_START'] + dist_to_mid
            return row['READ_START'] - tss_pos
        elif row['READ_STRAND'] == '-':
            tss_pos = row['TSS_END'] - dist_to_mid
            return tss_pos - row['READ_END'] 
        else:
            sys.exit('Encountered a read with no strand.')

    table['OFFSET_FROM_START'] = table.apply(offset_from_start, axis=1)

    # get three prime of cleavage sequences
    table['THREEPRIME_OF_CLEAVAGE'] = table.apply(lambda x: get_threeprime_sequence(x, genome_dict), axis=1)

    # reset index
    table = table.reset_index()
    # get the number of times it maps to a zero position
    zeros = table[table['OFFSET_FROM_START'] == 0].groupby('READ_NAME').size()
    unique_reads = table['READ_NAME'].unique()
    zeros = zeros.reindex(index=unique_reads, fill_value=0.)

    # get the zeros to total ratio
    zeros_to_total = zeros.astype(float) / num_times_map.astype(float)

    # set the index back to read name so we can fill in these values
    table = table.set_index('READ_NAME')
    
    # Figure out the ratio of reads mapping to the zero position of the tss
    # vs. nonzero for each read
    for vector in ('ZEROS', 'ZEROS_TO_TOTAL'):
        table[vector] = table[vector].fillna(eval(vector.lower()))

    # reset index
    table = table.reset_index()

    # Abundances per sequence (Series: sequence -> abundance)
    abundances = get_sequence_abundances(table)
    table = table.set_index('ORIGINAL_SEQUENCE')
    table['SEQUENCE_COUNT'] = table['SEQUENCE_COUNT'].fillna(abundances)
    table = table.reset_index()

    # reset columns if necessary (is it necessary?)
    table = table.reindex_axis(BEDTOOLS_INTERSECT_COLUMNS + COLUMNS_TO_ADD + \
                               GTF_KEYS, axis=1)
    # Sorting takes too long? 
    # table = table.sort_index(by='sequence_count', ascending=False)

    # Output the file
    table.to_csv(args.outfile, sep='\t', index=False, na_rep='NA')

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-f', '--fasta-file', type=str, nargs=1,
                        required=True, help=('hg19.fa; must be fasta file.'))
    parser.add_argument('-r', '--original-reads-file', type=str, nargs=1,
                        required=True, help=('$(INT_DIR)/%_before_mapping.fastq;'
                                             ' must be fastq file.'))
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args(argv)
    process_files(args)
    
if __name__ == '__main__':
    main(sys.argv[1:])
