import argparse
import pandas as pd
from pandas import DataFrame, Series
import sys
from mpltools import style
import matplotlib.pyplot as plt

def calculate_num_times_map(df):
    '''Adds a column, num_times_map, to the existing DataFrame.'''
    num_times_map = df.groupby('READ_NAME').size()
    df = df.set_index('READ_NAME')
    df['NUM_TIMES_MAP'] = num_times_map
    df['NUM_TIMES_MAP'][df['TSS_CHROM'] == '*'] = 0
    df = df.reset_index()
    return df

def convert_single_g(row):
    if row['3PTRIMMED'] == 'G':
        row['3PTRIMMED'] = None
        row['ORIGINAL_SEQUENCE'] += 'G'
    return row

def main(argv):
    style.use('ggplot')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=str, required=True)
    parser.add_argument('--stats-file', type=argparse.FileType('w'))
    parser.add_argument('--threeprime-trimmed', type=str)
    parser.add_argument('--three-prime-nucs', type=int, default=1)
    parser.add_argument('--convert-single-g', action='store_true')
    args = parser.parse_args(argv)
    df = pd.read_table(args.infile, index_col=None)
    if args.convert_single_g:
        df = df.apply(convert_single_g, axis=1)
    df = df[df['OFFSET_FROM_START'] == 0]
    if args.threeprime_trimmed is not None:
        if args.threeprime_trimmed == '':
            df = df[df['3PTRIMMED'].isnull()]
        else:
            df = df[df['3PTRIMMED'] == args.threeprime_trimmed]
    df = calculate_num_times_map(df)
    df = df.apply(lambda x: pd.Series(
        [x['ORIGINAL_SEQUENCE'][-1].upper() + \
        x['THREEPRIME_OF_CLEAVAGE'][0:args.three_prime_nucs].upper(),
        float(x['COUNT']) / float(x['NUM_TIMES_MAP'])]), axis=1)
    s = df.groupby(0).sum().sort(columns=1, ascending=False)
    s.plot(kind='bar')
    ax = plt.gca()
    ax.legend().set_visible(False)
    plt.savefig(args.outfile, format='eps')
    s = s.reset_index()
    total = s[1].sum()
    gs_pos_one = s[s[0].map(lambda x: x[1] == 'G')][1].sum()
    as_pos_zero = s[s[0].map(lambda x: x[0] == 'A')][1].sum()
    gs_exclusive = s[s[0].map(lambda x: x[1] == 'G' and x[0] != 'A')][1].sum()
    as_exclusive = s[s[0].map(lambda x: x[0] == 'G' and x[1] != 'G')][1].sum()
    gs_or_as_pos_zero = s[s[0].map(lambda x: x[0] == 'A' or x[0] == 'G')][1].sum()
    gs_both_exclusive = s[s[0].map(lambda x: 'G' in x and x[0] != 'A')][1].sum()
    gs_or_as_pos_zero_exclusive = s[s[0].map(lambda x: (x[0] == 'A' or x[0] == 'G') and x[1] != 'G')][1].sum()
    gs_pos_zero = s[s[0].map(lambda x: x[0] == 'G')][1].sum()
    if args.stats_file:
        for frac, string in [
                (gs_pos_one / total, 'Fraction Gs at position 1: %f\n'),
                (as_pos_zero / total, 'Fraction As at position 0: %f\n'),
                (gs_exclusive / total, 'Fraction Gs at position 1 with no As: %f\n'),
                (as_exclusive / total, 'Fraction As at position 0 with no Gs: %f\n'),
                (gs_or_as_pos_zero / total, 'Fraction As or Gs at pos 0: %f\n'),
                (gs_pos_zero / total, 'Fraction Gs at position 0: %f\n'),
                ((as_exclusive / (gs_both_exclusive + as_exclusive)), 'Fraction A over total: %f\n'),
                (total, 'Total: %f\n')
        ]:
            args.stats_file.write(string % frac)
    
if __name__ == '__main__':
    main(sys.argv[1:])
