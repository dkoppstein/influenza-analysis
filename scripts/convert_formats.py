from Bio import SeqIO
import argparse
import sys

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('--in-format', default='fastq')
    parser.add_argument('--out-format', default='fastq-illumina')
    args = parser.parse_args(argv)
    for line in SeqIO.parse(args.infile, format=args.in_format):
        args.outfile.write(line.format(args.out_format))

if __name__ == '__main__':
    main(sys.argv[1:])
