'''
Discard sequences greater than HIGH and/or shorter than LOW. 
'''

import Bio.SeqIO
import argparse
import sys

def process_file(args):
    for record in Bio.SeqIO.parse(args.infile, args.format):
        length = len(record.seq)
        if length >= args.low and length <= args.high:
            args.outfile.write(record.format(args.format))

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', default=sys.stdin,
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', default=sys.stdout,
                        type=argparse.FileType('w'))
    parser.add_argument('-g', '--high', type=int, default=sys.maxint)
    parser.add_argument('-l', '--low', type=int, default=-1)
    parser.add_argument('-f', '--format', type=str, default='fastq-illumina')
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
