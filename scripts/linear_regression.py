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


TIME_POINTS = [30, 45, 60, 90, 120, 240] # minutes

def heatmap(df, outfile, color_map='coolwarm', format='eps'):
    fig, ax = plt.subplots(figsize=(20,10))
    c = mcolors.ColorConverter().to_rgb
    cmap = plt.get_cmap(color_map)
    df_out = df.iloc[:, 1:-1] # grab everything but sequence and slope
    sys.stderr.write(str(df_out))
    hmap = ax.pcolor(df_out, cmap=cmap)
    ax.autoscale(tight=True)  # get rid of whitespace in margins of heatmap
    # ax.set_aspect('equal')  # ensure heatmap cells are square
    ax.xaxis.set_ticks_position('top')  # put column labels at the top
    ax.set_yticks(np.arange(len(df.index)) + 0.5)
    ax.set_yticklabels(df['SEQUENCE'], size=20)
    plt.xlabel('Time Course (minutes)', fontsize=36)
    ax.set_xticks(np.arange(len(df_out.columns)) + 0.5)
    ax.set_xticklabels(TIME_POINTS, size=24)

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "3%", pad="1%")
    cbar = plt.colorbar(hmap, cax=cax)
    cbar.set_label(r'Counts per million', fontsize=36)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(20)

    plt.savefig(outfile, format=format)
    
def process_df(final_df):
    # drop where any data is missing
    final_df.dropna(inplace=True)
    
    # collect the first coefficient of each row
    coeffs = final_df.drop('SEQUENCE', 1).apply(
        lambda x: np.polyfit(TIME_POINTS, x, 1)[0], axis=1)
    coeffs = DataFrame(coeffs)
    coeffs.columns = ['SLOPE']
    final_df = pd.concat([final_df, coeffs], axis=1)
    final_df = final_df.sort(columns=['SLOPE'], ascending=True)
    final_df.index = range(len(final_df)) # reindex
    return final_df
    
def process_file(args):
    final_df = None
    for time_point in TIME_POINTS:
        fn = '%s%s%s' % (args.infile_prefix, time_point, args.infile_suffix)
        if args.gzipped:
            fn = StringIO(Popen(['zcat', fn],
                        stdout=PIPE).communicate()[0])
        df = pd.read_table(fn, index_col=None)
        cname = '%s_count' % time_point
        df.columns = ['SEQUENCE', cname]
        if args.normalize:
            df[cname] /= float(df[cname].sum()) # counts per million
            df[cname] *= 1000000.
        if final_df is None:
            final_df = df
        else:
            final_df = pd.merge(final_df, df, on='SEQUENCE', how='outer', sort=False)
    if args.quantile_normalize:
        final_df = quantile_normalize(final_df,
            columns=[c for c in final_df.columns if c != 'SEQUENCE'])
    final_df = process_df(final_df)
    if args.first is not None:
        final_df = final_df[args.first:]
    if args.last is not None:
        final_df = final_df[:-args.last]
    heatmap(final_df, args.outfile, format=args.format)

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gzipped', action='store_true')
    parser.add_argument('-p', '--infile-prefix', type=str, required=True)
    parser.add_argument('-s', '--infile-suffix', type=str, required=True)
    parser.add_argument('--first', type=int, default=None)
    parser.add_argument('--last', type=int, default=None)
    parser.add_argument('--normalize', action='store_true')
    parser.add_argument('--quantile-normalize', action='store_true')
    parser.add_argument('-f', '--format', type=str, default='eps')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args(argv)
    return process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
