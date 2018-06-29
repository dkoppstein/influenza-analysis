
'''
Attempts to offset the intervals of a GTF file such that no Gs are at the front of each read. 
'''

import sys
from pybedtools import BedTool
import argparse
from Bio.Seq import Seq

def process_file(args):
    gtf = BedTool(args.infile[0])
    for iv in gtf:
        seq = BedTool.seq((iv.chrom, iv.start, iv.end), args.fasta_file[0])
        if iv.strand == '-':
            seq = str(Seq(seq).reverse_complement())
        num_gs = 0
        while seq.upper().startswith('G'):
            num_gs += 1
            seq = seq[1:]
        if iv.strand == '+':
            iv.start += num_gs
        elif iv.strand == '-':
            iv.end -= num_gs
        args.outfile.write(str(iv))

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=str, nargs=1,
                        required=True)
    parser.add_argument('-f', '--fasta-file', type=str, nargs=1,
                        required=True)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args(argv)
    process_file(args)

if __name__ =='__main__':
    main(sys.argv[1:])
