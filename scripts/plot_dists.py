from cStringIO import StringIO
import numpy as np
import pandas as pd
from pandas import Series, DataFrame
import sys
import argparse
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib import cm
from subprocess import Popen, PIPE
from flulib import quantile_normalize
from do_all_comparisons import DATASETS
import seaborn as sns

def findnth(haystack, needle, n):
    parts= haystack.split(needle, n+1)
    if len(parts)<=n+1:
        return -1
    return len(haystack)-len(parts[-1])-len(needle)

def get_basename(filename):
    if 'intermediate' in filename:
        filename = filename[len('intermediate')+1:]
    if filename.startswith('PELCHAT'):
        underscore = 3
    else:
        underscore = 1
    return filename[:findnth(filename, '_', underscore)]

def process_file(args):
    global DATASETS
    DATASETS = dict(DATASETS)
    final_df = None
    for f in args.infiles:
        name = f.name
        if args.gzipped:
            f = StringIO(Popen(['zcat', f.name],
                        stdout=PIPE).communicate()[0])
        df = pd.read_table(f, index_col=None)
        df.columns = ['SEQUENCE', name]
        if args.cutoff_reads:
            df = df[df[name] >= args.cutoff_reads]
        if args.normalize or args.normalize_to_one:
            df[name] /= float(df[name].sum()) # counts per million
            if not args.normalize_to_one:
                df[name] *= 1000000.
        if args.log_x:
            df[name] = df[name].map(np.log10)
        if final_df is None:
            final_df = df
        else:
            final_df = pd.merge(final_df, df, on='SEQUENCE',
                                how='outer', sort=False)
    
    cols = [col for col in final_df.columns if col != 'SEQUENCE']
    if args.plot_lengths:
        final_df['SEQUENCE'] = final_df['SEQUENCE'].map(len)
        final_df = final_df.groupby('SEQUENCE').sum()
    sns.set_style("darkgrid", {"axes.facecolor": "#E5E5E5", 'font.family': 'Arial'})
    sns.set_context('paper')
    with sns.color_palette('husl', 8, desat=0.6):
        fig, ax = plt.subplots(figsize=(10,10))
        ax.grid(False, which='minor')
        ax.grid(True, which='major')
        for col in cols:
            s = final_df[col].dropna()
            sns.distplot(s, ax=ax, bins=args.bins,
                         kde_kws={'shade': False,
                                  'label': DATASETS[get_basename(col)]},
                         hist_kws={'histtype': 'stepfilled'})
        plt.xlabel(args.xlabel)
        plt.ylabel(args.ylabel)
        plt.savefig(args.outfile, format=args.format, dpi=args.dpi)

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gzipped', action='store_true')
    parser.add_argument('-i', '--infiles', type=argparse.FileType('r'),
                        required=True, nargs='+')
    parser.add_argument('--normalize', action='store_true')
    parser.add_argument('--normalize-to-one', action='store_true')
    parser.add_argument('--bins', default=20, type=int)
    parser.add_argument('--quantile-normalize', action='store_true')
    parser.add_argument('--plot-lengths', action='store_true')
    parser.add_argument('-f', '--format', type=str, default='png')
    parser.add_argument('--xlabel', type=str, default='')
    parser.add_argument('--cutoff-reads', default=None, type=int)
    parser.add_argument('--ylabel', type=str, default='')
    parser.add_argument('-o', '--outfile', type=str)
    parser.add_argument('--dpi', default=300, type=int)
    parser.add_argument('--log-x', action='store_true')
    args = parser.parse_args(argv)
    return process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
