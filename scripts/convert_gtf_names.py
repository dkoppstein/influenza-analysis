import argparse
import pandas as pd
import sys
from preprocess_bedtools_output import deparse_gtf_metadata

'Converts the NAME column of a BED file to the GTF name specified in the last column (i.e. transcript_id)'

def process_file(args):
    df = pd.read_table(args.infile, index_col=None, header=None)
    metaseries = df[df.columns[-1]].map(deparse_gtf_metadata)
    df[3] = metaseries.map(lambda x: x[args.to.upper()])
    df.to_csv(args.outfile, index=False, header=None, sep='\t')

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('--to', type=str, default='transcript_id')
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
