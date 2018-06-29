
import argparse
import sys
from Bio import SeqIO
import pandas as pd

def process_files(args):
    cdf = pd.read_table(args.cufflinks_file, index_col=None)
    cdf = cdf.sort(columns=['FPKM'], ascending=False)
    cdf = cdf.head(args.num_genes)
    genes = set(cdf['tracking_id'])
    for record in SeqIO.parse(args.fasta_file, 'fasta'):
        if record.id in genes:
            args.outfile.write(record.format('fasta'))

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta-file', type=argparse.FileType('r'))
    parser.add_argument('-c', '--cufflinks-file', type=argparse.FileType('r'),
                        help='genes.fpkm_tracking')
    parser.add_argument('--num-genes', type=int, default=5000)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args(argv)
    process_files(args)

if __name__ == '__main__':
    main(sys.argv[1:])
