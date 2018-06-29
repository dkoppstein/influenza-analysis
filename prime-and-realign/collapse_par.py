#!/usr/bin/env python
# collapse_par.py

__author__ = 'David Koppstein'

'''
Takes a fastq file of collapsed sequences which contain essentially only 
the cellular mRNA 5' end, as well as a 5' end of an influenza transcript. 

Outputs a fastq file to stdout with all possible prime-and-realign sequences
removed. NB: if the subsequence is not present, outputs zero. 

Prints statistics to --stats-file. 
'''

import argparse
import sys
from Bio import SeqIO
from collections import defaultdict
import re
import pandas as pd

def add_metadata_to_header(record, ann_dict, delim='|', eq='='):
    annotation = ''
    for k, v in ann_dict.items():
        annotation += '{}{}{}{}'.format(delim, k.upper(), eq, v)
    record.description += annotation
    record.id = ''


def process_file(args):
    # quick access to some helpful args
    infile, sequence = args.infile, args.sequence.upper()
    
    # dictionary: string -> int
    # d.key was removed d.item times
    d = defaultdict(int)
    max_length = args.max_length

    # do the analysis
    for record in SeqIO.parse(args.infile, args.format):
        final_seq_trimmed = ''
        for i in range(args.num_trims): # number of times to trim pars
            length = args.max_length
            while length >= args.min_length:
                seq = sequence[:length]
                if str(record.seq).upper().endswith(seq):
                    record = record[:-length]
                    final_seq_trimmed = '%s%s' % (seq, final_seq_trimmed)
                    break
                else:
                    length -= 1
        if args.header:
            add_metadata_to_header(record, {'3PTRIMMED': final_seq_trimmed})
        d[final_seq_trimmed] += 1
        args.outfile.write(record.format(args.format))
                
    # Write statistics to stderr
    if args.stats_file:
        pd.Series(d).fillna(0).to_csv(args.stats_file, sep='\t')


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', default=sys.stdin,
                        type=argparse.FileType('r'),
                        help='A file of BioPython SeqRecords file with only cellular 5 prime ends')
    parser.add_argument('-f', '--format', default='fastq')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-s', '--sequence', type=str, help=
                        'Sequence of the 5\' end of influenza',
                        required=True)
    parser.add_argument('-m', '--max-length', type=int, default=10,
                        help='The max length of possible prime-and-realign events')
    parser.add_argument('-n', '--min-length', type=int, default=0, 
                        help='The min length of possible prime-and-realign events')
    parser.add_argument('-d', '--header', action='store_true', 
                        default=False, help='If header, adds " ' + \
                        '|3PTRIMMED=SEQUENCE| to the header.')
    parser.add_argument('-t', '--stats-file', type=argparse.FileType('w'),
                        required=False)
    parser.add_argument('--num-trims', type=int, help='Number of times '
                        'to trim a sequence.', required=True)
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
