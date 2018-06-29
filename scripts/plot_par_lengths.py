import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from mpltools import style
import seaborn as sns
from matplotlib.pyplot import figaspect
from matplotlib.figure import Figure

# Takes a set of par_stats files and plots a histogram of the lengths

WIDTH = 0.35

# convenience function for sizing
def cm2inch(*tupl):
    inch = 2.54
    if type(tupl[0]) == tuple:
        return tuple(i/inch for i in tupl[0])
    else:
        return tuple(i/inch for i in tupl)

def row_apply(row):
    if type(row['SEQUENCE']) is str:
        length = len(row['SEQUENCE'])
    else:
        length = 0
    return pd.Series([length, row['OCCURRENCES']])

def process_files(args):
    sns.set_context("paper", rc={"lines.linewidth": 0.4})
    sns.set_style("darkgrid", {"axes.facecolor": "#E5E5E5", 'font.family': 'Arial'})

    w, h = figaspect(1.)
    fig = Figure(figsize=(w,h))
    ax = plt.gca()
#    fig.set_size_inches(*cm2inch(8.2, 8.2))

    for file_list in [args.GCAAAAGCAG_files, args.GCGAAAGCAG_files]:
        final_series = None
        for handle in file_list:
            df = pd.read_table(handle, index_col=None)
            df2 = df.apply(row_apply, axis=1)
            df2.columns = ['LENGTH', 'COUNT']
            series = df2.groupby('LENGTH').sum()
            series = series / series.sum()
            series = series[args.min_num:args.max_num]
            if final_series is None:
                final_series = series
            else:
                final_series = pd.concat([final_series, series], axis=1)
        stdevs = final_series.apply(lambda x: np.std(x), axis=1)
        if file_list is args.GCAAAAGCAG_files:
            rects1 = ax.bar(series.index, series['COUNT'], WIDTH, 
                            color='r', yerr=np.array(stdevs),
                            ecolor='black')
        else:
            rects2 = ax.bar(series.index + WIDTH, series['COUNT'], WIDTH,
                            color='b', yerr=np.array(stdevs),
                            ecolor='black')

    ax.legend((rects1[0], rects2[0]), ('HA, NA, NP, NS1', 'MP, PA, PB1, PB2'), fontsize=16)
    ax.set_xticks(series.index + WIDTH)
    ax.set_xticklabels(series.index)
    ax.set_xlabel('Nucleotides Matching Viral Genome', fontsize=16)
    ax.set_ylabel('Frequency', fontsize=16)
    
    for axis in (ax.xaxis, ax.yaxis):
        for tick in axis.get_major_ticks():
            tick.label.set_fontsize(16)

    plt.axes().set_aspect(4 / 0.35)            
    plt.tight_layout()
    plt.savefig(args.outfile, format=args.format)

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--GCAAAAGCAG-files', nargs='+',
                        type=argparse.FileType('r'))
    parser.add_argument('--GCGAAAGCAG-files', nargs='+',
                        type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', type=str, required=True)
    parser.add_argument('--min-num', type=int, default=2)
    parser.add_argument('--max-num', type=int, default=10)
    parser.add_argument('-f', '--format', type=str, default='eps')
    args = parser.parse_args(argv)
    process_files(args)

if __name__ == '__main__':
    main(sys.argv[1:])
