#!/usr/bin/env python
# collapse_par.py

__author__ = 'David Koppstein'

'''
Takes a fastq file of collapsed sequences which contain essentially only 
the cellular mRNA 5' end, as well as a 5' end of an influenza transcript. 

Outputs a fastq file to stdout with all possible prime-and-realign sequences
removed. NB: if the subsequence is not present, outputs zero. 

Prints statistics to stderr. 
'''

FORMAT = 'fastq-illumina'
COUNT_PAT = 'count=(\d+)'

import flulib
import argparse
import sys
import Bio.SeqIO
from collections import defaultdict
import re

def get_field(record_id, prefix):
    pat = prefix + '=([-]?\d+)'
    return re.search(pat, record_id).groups()[0]

def get_count(record_id):
    'Return an int of the count in the record header'
    pat = 'count=(\d+)'
    return int(re.search(pat, record_id).groups()[0])

def process_file(args):
    # quick access to some helpful args
    infile, sequence = args.infile, args.sequence.upper()
    
    # dictionary: string -> int
    # d.key was removed d.item times
    d = defaultdict(int)
    max_length = args.max_length

    # prepopulate the dict with zeros
    # don't use a default dict because we actually want the zeros
    while max_length > 0:
        d[sequence[0:max_length]] = 0
        max_length -= 1

    # do the analysis
    for record in Bio.SeqIO.parse(args.infile, FORMAT):
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
            flulib.add_metadata_to_header(record,
                                        {'3PTRIMMED': final_seq_trimmed})
        if args.count:
            d[final_seq_trimmed] += get_count(record.id)
        else: 
            d[final_seq_trimmed] += 1
        args.outfile.write(record.format(FORMAT))
                
    # Write statistics to stderr
    if args.stats_file:
        args.stats_file.write('SEQUENCE\tOCCURRENCES\n')
        for k, v in d.iteritems():
            args.stats_file.write(k + '\t' + str(v) + '\n')

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', default=sys.stdin,
                        type=argparse.FileType('r'),
                        help='A Solexa 1.3-1.7 fastq file with only '
                        'cellular 5 prime ends')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-s', '--sequence', type=str, help=
                        'Sequence of the 5\' end of influenza',
                        required=True)
    parser.add_argument('-c', '--count', action='store_true', 
                        default=False, help='If --count, use the header ' + \
                        'count=N in the FASTQ file to determine count.')
    parser.add_argument('-m', '--max-length', type=int, default=10,
                        help='The max length of possible prime-and-realign events')
    parser.add_argument('-n', '--min-length', type=int, default=0, 
                        help='The min length of possible prime-and-realign events')
    parser.add_argument('-d', '--header', action='store_true', 
                        default=False, help='If header, adds " ' + \
                        '|3ptrimmed=SEQUENCE| to the header.')
    parser.add_argument('-f', '--stats-file', type=argparse.FileType('w'),
                        required=False)
    parser.add_argument('--num-trims', type=int, help='Number of times '
                        'to trim a sequence.', required=True)
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
