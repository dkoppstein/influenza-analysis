import pandas as pd
import argparse
import sys

def process_seq(s, lower=None, upper=None):
    if lower:
        if len(s) < lower:
            return ''
    else:
        lower = 0
    if (upper and upper > len(s)) or not upper:
        higher = len(s)
    return s[lower:upper]

def process(args):
    df = pd.read_table(args.infile, index_col=None)
    if type(df[df.columns[0]][0]) is str:
        seq_col, weight_col = df.columns[0], df.columns[1]
    else:
        seq_col, weight_col = df.columns[1], df.columns[0]
    df[seq_col] = df[seq_col].map(
        lambda x: process_seq(x, lower=args.lower, upper=args.upper))
    df.to_csv(args.outfile, sep='\t', index=False)
    

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-l', '--lower', type=int, default=None)
    parser.add_argument('-u', '--upper', type=int, default=None)
    args = parser.parse_args(argv)
    process(args)

if __name__ == '__main__':
    main(sys.argv[1:])
