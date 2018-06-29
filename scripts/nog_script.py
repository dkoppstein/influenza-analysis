# quick script to calculate number of reads that start with a G

import sys
from Bio import SeqIO
import argparse

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-f', '--format', default='fastq-illumina',
                        type=str)
    args = parser.parse_args(argv)
    g = 0
    nog = 0
    for record in SeqIO.parse(args.infile, args.format):
        if str(record.seq).upper().endswith('G'):
            g += 1
        else:
            nog += 1
    args.outfile.write('Gs\tNo_Gs\n')
    args.outfile.write('%s\t%s\n' % (g, nog))

if __name__ == '__main__':
    main(sys.argv[1:])
