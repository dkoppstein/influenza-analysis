from __future__ import absolute_import

import HTSeq
import sys
import pandas as pd
import flulib
from Bio import SeqIO
import argparse
import numpy as np
import time
if pd.__version__ < '0.13.1':
    raise ImportError('Requires pandas 0.13.1.')

from deparse_bedtools_output import (
    BEDTOOLS_INTERSECT_COLUMNS, COLUMNS_TO_ADD, GTF_KEYS
    )

'''
This file takes bedtools output and creates a single table with all
possible information about a set of reads.
'''
    
def get_sequence_abundances(table):
    '''
    Returns the abundance of each sequence, in counts, as a pandas Series
    sequence -> count
    '''
    no_dups = table.drop_duplicates(cols='READ_NAME', take_last=True)
    return no_dups['COUNT'].astype(int).groupby(no_dups['ORIGINAL_SEQUENCE']).sum()
    
def process_files(args):
    # Add reads that didn't map back to the table (bamToBed filters them out!!)
    reader = HTSeq.BAM_Reader(args.bam_file[0])
    unmapped_reads = []
    for aln in reader:
        if not aln.aligned: # select only unmapped reads
            read_name_dict = flulib.headerdict_from_string(aln.read.name)
            read_name_dict = {k.upper() : v for k, v in read_name_dict.iteritems()}
            temp_d = {
                'READ_NAME': aln.read.name,
                'ORIGINAL_SEQUENCE': aln.read.seq,
                'ORIGINAL_SEQUENCE_LENGTH': len(aln.read.seq),
                'NUM_TIMES_MAP': 0,
                'TSS_CHROM': '.',
                'ZEROS': 0,
            }
            temp_d.update(read_name_dict)
            unmapped_reads.append(temp_d)
            for k in BEDTOOLS_INTERSECT_COLUMNS + COLUMNS_TO_ADD + GTF_KEYS:
                if k not in temp_d:
                    temp_d[k] = None
    
    # Create a new DataFrame with unmapped reads
    # and concatenate it with the table
    table = pd.DataFrame(unmapped_reads)

    # Abundances per sequence (Series: sequence -> abundance)
    abundances = get_sequence_abundances(table)
    table = table.set_index('ORIGINAL_SEQUENCE')
    table['SEQUENCE_COUNT'] = table['SEQUENCE_COUNT'].fillna(abundances)
    table = table.reset_index()

    # Reset columns if necessary
    table = table[BEDTOOLS_INTERSECT_COLUMNS + COLUMNS_TO_ADD + GTF_KEYS]
    # Sort by sequence count for easy visualization
    # Sorting takes too long? 
    # table = table.sort_index(by='sequence_count', ascending=False)
    
    # Output the file
    table.to_csv(args.outfile, sep='\t', index=False, na_rep='NA')

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-b', '--bam-file', type=str, nargs=1)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args(argv)
    process_files(args)
    
if __name__ == '__main__':
    main(sys.argv[1:])
