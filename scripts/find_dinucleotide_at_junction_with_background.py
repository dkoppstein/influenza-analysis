import argparse
import pandas as pd
from pandas import DataFrame, Series
import sys
from scipy import log2
from mpltools import style
import matplotlib.pyplot as plt
from cStringIO import StringIO
from subprocess import Popen, PIPE
from collections import defaultdict
from itertools import product
import numpy as np
import seaborn as sns

def calculate_num_times_map(df):
    '''Adds a column, num_times_map, to the existing DataFrame.'''
    num_times_map = df.groupby('READ_NAME').size()
    df = df.set_index('READ_NAME')
    df['NUM_TIMES_MAP'] = num_times_map
    df['NUM_TIMES_MAP'][df['TSS_CHROM'] == '*'] = 0
    df = df.reset_index()
    return df

def calculate_background_dinucleotides(df):
    if type(df[df.columns[0]][0]) is str:
        df.columns = ['SEQUENCE', 'WEIGHT']
    else:
        df.columns = ['WEIGHT', 'SEQUENCE']
    d = defaultdict(int)
    for idx, row in df.iterrows():
        s = row['SEQUENCE']
        if len(s) > 1:
            for i in range(len(s) - 1):
                d[s[i:i+2].upper()] += row['WEIGHT']
    return pd.Series(d)

def convert_single_g(row):
    if row['3PTRIMMED'] == 'G':
        row['3PTRIMMED'] = None
        row['ORIGINAL_SEQUENCE'] += 'G'
    return row

def main(argv):
    sns.set_style("darkgrid", {"axes.facecolor": "#E5E5E5", 'font.family': 'Arial'})
    sns.set_context('paper')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=str, required=True)
    parser.add_argument('--background-file', type=str,
                        default=None)
    parser.add_argument('--gzipped-background-file', action='store_true')
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
        float(x['WEIGHT'])]), axis=1)
    series = df.groupby(0).sum()
    series = pd.Series(series[1], index=series.index)
    nucs = ['A', 'C', 'G', 'T']
    all_dinucs = [x[0] + x[1] for x in product(nucs, nucs)]
    final_series = pd.Series([None]*len(all_dinucs),
                             index=all_dinucs).fillna(series).fillna(0)
    final_series = final_series.astype(float) / final_series.sum()
    if args.background_file is not None:
        if args.gzipped_background_file:
            args.background_file = StringIO(
                Popen(['zcat', args.background_file],
                    stdout=PIPE).communicate()[0])
        background = pd.read_table(args.background_file, index_col=None)
        background = calculate_background_dinucleotides(background)
        final_background = pd.Series(
            [None]*len(all_dinucs),
            index=all_dinucs).fillna(background).fillna(0)
        final_background = final_background.astype(float) / final_background.sum()
        final_series = final_series.astype(float) / \
          final_background.astype(float)
    final_series = final_series.replace([np.inf, -np.inf], np.nan)
    final_series = final_series.fillna(0.)
    final_series = final_series.map(log2)
    final_series.sort(ascending=False)
    final_series.plot(kind='bar')
    ax = plt.gca()
    ax.legend().set_visible(False)
    plt.ylabel('log$_2$(enrichment)')
    plt.ylim(-2.0, 2.0)
    plt.savefig(args.outfile, format='eps')
    final_series = final_series.reset_index()
    total = final_series.sum()
    gs_pos_one = final_series[final_series.index.map(lambda x: x[1] == 'G')][1].sum()
    as_pos_zero = final_series[final_series[0].map(lambda x: x[0] == 'A')][1].sum()
    gs_exclusive = final_series[final_series[0].map(lambda x: x[1] == 'G' and x[0] != 'A')][1].sum()
    as_exclusive = final_series[final_series[0].map(lambda x: x[0] == 'G' and x[1] != 'G')][1].sum()
    gs_or_as_pos_zero = final_series[final_series[0].map(lambda x: x[0] == 'A' or x[0] == 'G')][1].sum()
    gs_both_exclusive = final_series[final_series[0].map(lambda x: 'G' in x and x[0] != 'A')][1].sum()
    gs_or_as_pos_zero_exclusive = final_series[final_series[0].map(lambda x: (x[0] == 'A' or x[0] == 'G') and x[1] != 'G')][1].sum()
    gs_pos_zero = final_series[final_series[0].map(lambda x: x[0] == 'G')][1].sum()
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
