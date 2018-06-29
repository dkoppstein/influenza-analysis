import sys
import pandas as pd
from Bio import SeqIO
import argparse
import flulib
from collections import defaultdict

def process_file(args):
    lst = []
    for record in SeqIO.parse(args.infile, args.format):
        if args.use_header:
            d = flulib.headerdict_from_string(record.id)
            d = {k.upper(): v for k, v in d.iteritems()}
        else:
            d = {}
        d['SEQUENCE'] = str(record.seq).upper()
        if args.include_read_name:
            d['READ_NAME'] = record.id
        lst.append(d)
    df = pd.DataFrame(lst)
    df.columns = [h.upper() for h in df.columns]
    if 'COUNT' not in df.columns:
        df = pd.DataFrame(df.groupby(list(df.columns)).size()).reset_index()
        df.rename(columns={0: 'COUNT'}, inplace=True)
    if not type(df['COUNT'][0]) is int:
        df['COUNT'] = df['COUNT'].astype(int)
    if args.summarize:
        df = df['COUNT'].groupby(df['SEQUENCE']).sum().reset_index()
        df.rename(columns={0: 'COUNT'}, inplace=True)
    if args.sort:
        df = df.sort(columns=['COUNT'], ascending=False)
    df.columns = [s.upper() for s in df.columns]
    df.to_csv(args.outfile, sep='\t', index=False)
    return df

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-f', '--format', type=str, default='fastq-illumina')
    parser.add_argument('-s', '--sort', action='store_true')
    parser.add_argument('--summarize', action='store_true')
    parser.add_argument('--use-header', action='store_true')
    parser.add_argument('--include-read-name', action='store_true')
    args = parser.parse_args(argv)
    return process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
