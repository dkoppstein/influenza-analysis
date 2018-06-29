from cStringIO import StringIO
from subprocess import Popen, PIPE
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import numpy as np
import Bio.SeqIO
import sys
from collections import defaultdict
from mpltools import style
import seaborn as sns
from do_all_comparisons import DATASETS
DATASETS = dict(DATASETS)
from plot_dists import get_basename

def process_files(args):
    sns.set_context("paper", rc={"lines.linewidth": 0.4})
    sns.set_style("darkgrid", {"axes.facecolor": "#E5E5E5", 'font.family': 'Arial'})
    sns.set_palette('hls', 8)

    final_df = None
    fig = plt.figure()
    ax = plt.gca()
    for f in args.infiles:
        name = f.name
        if args.gzipped:
            f = StringIO(Popen(['zcat', name],
                        stdout=PIPE).communicate()[0])
        df = pd.read_table(f, index_col=None)
        df.columns = ['SEQUENCE', get_basename(name)]
        if final_df is None:
            final_df = df
        else:
            final_df = pd.merge(final_df, df, on='SEQUENCE',
                                how='outer', sort=False)
    cols = [col for col in final_df.columns if col != 'SEQUENCE']
    final_df['SEQUENCE'] = final_df['SEQUENCE'].map(len)
    final_df = final_df.groupby('SEQUENCE').sum()
    for col in cols:
        final_df[col] /= final_df[col].sum()
    for col in cols:
        s = final_df[col].dropna()
        plt.plot([0] + list(final_df.index), [0] + list(final_df[col]), label=DATASETS[col])
    plt.xlabel(args.xlabel, fontsize=16)
    plt.ylabel(args.ylabel, fontsize=16)
    plt.xticks(np.arange(0, args.max_x, 1.0))
    plt.legend(loc='best', fontsize=16)
    for axis in (ax.xaxis, ax.yaxis):
        for tick in axis.get_major_ticks():
            tick.label.set_fontsize(16)

    plt.xlim(0, args.max_x)
    plt.ylim(0, 0.35)
    ax.set_aspect(20 * 4 / 5.25 / 0.35) # empirically derived i hate matplotlib
    plt.tight_layout()
    plt.savefig(args.outfile, format=args.output_format)

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=\
                                     argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infiles', nargs='+',
                        type=argparse.FileType('r'), help='Input files.')
    parser.add_argument('--gzipped', action='store_true')
    parser.add_argument('-f', '--format', nargs='?', type=str,
                        default='fastq-illumina', 
                        choices=['fastq-illumina', 'fastq-solexa', 'fastq',
                                 'table', 'fasta'],
                        help='Format of input file.')
    parser.add_argument('-o', '--outfile', nargs='?', 
                        type=str, required=True)
    parser.add_argument('--xlabel', type=str, default='Length')
    parser.add_argument('--ylabel', type=str, default='Frequency')
    parser.add_argument('--max-x', type=int, default=20)
    parser.add_argument('-u', '--output-format', type=str, default='eps',
                        nargs='?')
    parser.add_argument('-t', '--title', type=str, default='', 
                        nargs='?')
    args = parser.parse_args(argv)
    process_files(args)

if __name__ == '__main__':
    main(sys.argv[1:])
