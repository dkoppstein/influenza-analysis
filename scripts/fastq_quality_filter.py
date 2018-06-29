import Bio.SeqIO
import argparse
import numpy as np
import sys

'''
Filters a FASTQ file according to the following criteria:
'''

def process_file(args):
    for record in Bio.SeqIO.parse(args.infile, args.format):
        if np.mean(record.letter_annotations.values()[0]) >= args.min_quality \
          and not (args.no_ns and 'N' in str(record.seq)):
            args.outfile.write(record.format(args.format))

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'), 
                        default=sys.stdin, help='Default is stdin.')
    parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help='Default is stdout.')
    parser.add_argument('--min-quality', type=int, default=0,
                        help='Minimum average sequence quality of the read.')
    parser.add_argument('--no-ns', action='store_true', 
                        help='If specified, will remove all reads containing Ns.')
    parser.add_argument('-f', '--format', default='fastq-illumina', nargs='?')
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
