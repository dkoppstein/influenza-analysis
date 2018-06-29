import pandas as pd
import sys
import argparse
from subprocess import Popen, PIPE
import matplotlib
import matplotlib.pyplot as plt
import prettyplotlib as ppl
from scipy.stats import spearmanr
import numpy as np
from cStringIO import StringIO
import seaborn as sns
from flulib import quantile_normalize
from matplotlib.markers import MarkerStyle

'Compares cufflinks-counted RNA-Seq files and subsets of data. '

SEQUENCE_COLORMAP = {
    'ATCGCTTCTCG': '#2078B4', # U2
    'ATCGCTTCTC': '#2078B4', # U2
    'ATACTTACCT': '#34A048', # U1
    'ATACTTACCTG': '#34A048', # # U1
    'AAGACTATATTT': '#CA4B9B', # U3
    'AAGACTATATTTT': '#CA4B9B', # U3
    'AAGACTATACTTTCA': '#CA4B9B', #U3
    'ATCGTCAGGT': '#708455', # U8
    'ATCGTCAGGTG': '#708455', # U8
    'ATCCTTTTGTA': '#F59597', # U13
    'AGCTTTGCGCA': '#E21F26', # U4
    'AGCTTTGC': '#E21F26', # U4
    'ATACTCTGGTTT': '#F79B2E', # U5
    'ATACTCTGGTTTCT': '#F79B2E', # U5
    'ATACTCTGGTTTCTTCA': '#F79B2E', # U5
    'AGTGTTACA': '#D75927', # U7
    'AAAAAGGGCTTCT': '#9378A3', # U11
    'AAAAAAGGGCTTCT': '#9378A3', # U11
    'ATGCCTTAAACTTA': '#5E479D', # U12
    'ATGCCTTAAACTTAT': '#5E479D' # U12
}
    
U_SNRNA_DATAFRAME = pd.DataFrame([
    ['ATCGCTTCTCG', '#2078B4', 'U2', 'o', '#2078B4'],
    ['ATCGCTTCTC', '#2078B4', 'U2', 'v', '#2078B4'],
    ['ATACTTACCT', '#34A048', 'U1', '^', '#34A048'],
    ['ATACTTACCTG', '#34A048', 'U1', '<', '#34A048'],
    ['AAGACTATATTT', '#CA4B9B', 'U3', '>', '#CA4B9B'],
    ['AAGACTATATTTT', '#CA4B9B', 'U3', 's', '#CA4B9B'],
    ['AAGACTATACTTTCA', '#CA4B9B', 'U3', 'D', '#CA4B9B'],
    ['ATCGTCAGGT', '#708455', 'U8', 'd', '#708455'],
    ['ATCGTCAGGTG', '#708455', 'U8', 'o', 'none'],
    ['ATCCTTTTGTA', '#F59597', 'U13', 'v', 'none'],
    ['AGCTTTGCGCA', '#E21F26', 'U4', '^', 'none'],
    ['AGCTTTGCGCAGTG', '#E21F26', 'U4', '<', 'none'],
    ['AGCTTTGCGCAGT', '#E21F26', 'U4', '>', 'none'],
    ['AGCTTTGC', '#E21F26', 'U4', 's', 'none'],
    ['ATACTCTGGTTT', '#F79B2E', 'U5', 'o', 'none'],
    ['ATACTCTGGTTTCT', '#F79B2E', 'U5', 'd', 'none'],
    ['ATACTCTGGTTTCTCTTCA', '#F79B2E', 'U5', 'o', '#F79B2E'],
    ['AGTGTTACA', '#6DCDE9', 'U7', 'D', 'none'],
    ['AAAAAGGGCTTCT', '#9378A3', 'U11', '^', '#9378A3'],
#    ['AAAAAAGGGCTTCT', '#9378A3', 'U11', '<', '#9378A3'],
    ['ATGCCTTAAACTTA', '#5E479D', 'U12', '>', '#5E479D'],
    ['ATGCCTTAAACTTAT', '#5E479D', 'U12', 's', '#5E479D']
], columns=['SEQUENCE', 'COLOR', 'NAME', 'SHAPE', 'FILL'])

def has_prime_and_realigned(original_seq, target_seq, par_sequence):
    return original_seq.startswith(target_seq) and \
        par_sequence.startswith(original_seq[len(target_seq):])
        
def process_files(args):
    sns.set_style("darkgrid", {"axes.facecolor": "#E5E5E5", 'font.family': 'Arial'})
    sns.set_context('paper')
    if args.no_header:
        names = ['SEQUENCE', 'COUNT']
    else:
        names = None
    if args.include_par:
        assert args.include_par == 'GCAAAAGCAG' or args.include_par == 'GCGAAAGCAG'
    if args.gzipped:
        args.file_a, args.file_b = (StringIO(Popen(['zcat', f.name],
                                          stdout=PIPE).communicate()[0]) \
                                          for f in [args.file_a, args.file_b])
    file_a = pd.read_table(args.file_a, index_col=None, names=names)
    file_b = pd.read_table(args.file_b, index_col=None, names=names)
    file_list = [file_a, file_b]
    for idx, f in enumerate(file_list):
        for col in (args.x_col, args.y_col):
            if args.cpm:
                f['CPM'] = f[col] / float(f[col].sum())
                f['CPM']= f['CPM'] * 1000000.
            if args.greater_than_zero:
                file_list[idx] = file_list[idx][file_list[idx][col] > 0]
            if args.cutoff_reads:
                file_list[idx] = file_list[idx][file_list[idx][col] >= args.cutoff_reads]
            if args.cutoff_cpm:
                file_list[idx] = file_list[idx][file_list[idx]['CPM'] >= args.cutoff_cpm]
    n_x = file_a[args.x_col].astype(int).sum()
    n_y = file_b[args.y_col].astype(int).sum()
    u_x, u_y = (len(df) for df in (file_a, file_b))

    merged_df = pd.merge(file_list[0], file_list[1], on=args.on, how=args.how).fillna(0.)
    if args.x_col == args.y_col:
        args.x_col += '_x'
        args.y_col += '_y'
    if args.only_starts_with is not None:
        for word in args.only_starts_with:
            merged_df = merged_df[merged_df[args.on].map(
                lambda x: x.startswith(word))]
    n_ux = merged_df[args.x_col].astype(int).sum()
    n_uy = merged_df[args.y_col].astype(int).sum()
    for col in (args.x_col, args.y_col):
        merged_df = merged_df.replace([np.inf, -np.inf], np.nan).dropna(subset=[col], how="all")
        merged_df[col] = merged_df[col].map(lambda x: np.log10(x))
        merged_df = merged_df.replace([np.inf, -np.inf], np.nan).dropna(subset=[col], how="all")
    if args.restrict_lengths is not None:
        filt = merged_df['SEQUENCE'].apply(
            lambda x: len(x) >= args.restrict_lengths[0] and \
            len(x) <= args.restrict_lengths[1])
        original_df = merged_df.copy()
        merged_df = merged_df[filt]
    if args.quantile_normalize:
        merged_df = quantile_normalize(
            merged_df, columns=[args.x_col, args.y_col])
    if not args.add_histograms:
        fig = plt.figure()
        ax = plt.gca()
        plt.scatter(merged_df[args.x_col], merged_df[args.y_col],
                    edgecolor=sns.desaturate('black', 0.75), facecolors='none',
                    s=0.5)
#        ax.set_yscale('log')
#        ax.set_xscale('log')
    else:
        g = sns.JointGrid(args.x_col, args.y_col, merged_df)
        if args.bandwidth:
            kde_kws = {'bw': args.bandwidth}
        else:
            kde_kws = None
        g.plot_marginals(sns.distplot, hist=True, kde=True,
                         color='blue', bins=50, kde_kws=kde_kws)
        g.plot_joint(plt.scatter, edgecolor=sns.desaturate('black', 0.75), facecolors='none',
                     s=0.5)
        ax = g.ax_joint
#        g.ax_marg_x.set_xscale('log')
#        g.ax_marg_x.set_yscale('log')
#        g.ax_marg_y.set_yscale('log')
#        g.ax_marg_y.set_xscale('log')
        g.ax_marg_x.grid(False, which='both')
        g.ax_marg_y.grid(False, which='both')
        ax = g.ax_joint
        ax.grid(False, which='minor')
    ax.grid(True, which='major')
#    ax.axis('equal')
    plt.xlabel(args.xlabel, fontsize=16)
    plt.title(args.title)
    plt.ylabel(args.ylabel, fontsize=16)
    plt.xlim([0.5, np.log10(args.max_val)])
    plt.ylim([0.5, np.log10(args.max_val)])
    for axis in (ax.xaxis, ax.yaxis):
        for tick in axis.get_major_ticks():
            tick.label.set_fontsize(16)
    if args.snRNAs is not None:
        if 'all' in args.snRNAs:
            args.snRNAs = list(U_SNRNA_DATAFRAME['SEQUENCE'])
            args.snRNAs.sort()
        for sequence in args.snRNAs:
            row = U_SNRNA_DATAFRAME[U_SNRNA_DATAFRAME.apply(
                lambda x: x['SEQUENCE'] == sequence, axis=1)]
            facecolor = row['FILL'][row.index[0]]
            edgecolor = row['COLOR'][row.index[0]]
            marker = row['SHAPE'][row.index[0]]
            snrnas = [None, None]
            if args.include_par:
                f = original_df['SEQUENCE'].map(
                    lambda x: has_prime_and_realigned(x, sequence, args.include_par))
                sub_df = original_df[f]
            if sequence == 'AGCTTTGCGCA':
                for i, label in enumerate([args.xlabel, args.ylabel]):
                    trimmed = False
                    for gene in ('HA', 'NA', 'NP', 'NS1'):
                        if gene in label:
                            trimmed = True
                            break
                    if trimmed:
                        snrnas[i] = file_list[i][file_list[i]['SEQUENCE'] == 'AGCTTT'][args.x_col[:-2]].map(lambda x: np.log10(x))
                    else:
                        snrnas[i] = file_list[i][file_list[i]['SEQUENCE'] == sequence][args.y_col[:-2]].map(lambda x: np.log10(x))
            elif sequence == 'ATACTCTGGTTTCTCTTCA' and 'NS1' in args.xlabel and \
            'PB2' in args.ylabel and not ('Pelchat' in args.xlabel):
                    # find numbers in vfiltered_summary.txt.gz file
                    snrnas[0] = pd.Series(np.log10(130.))
                    snrnas[1] = pd.Series(np.log10(12.))
            else:
                sub_df = merged_df[merged_df['SEQUENCE'] == sequence]
                snrnas[0] = sub_df[args.x_col]
                snrnas[1] = sub_df[args.y_col]
            if len(snrnas[0]) > 0 and len(snrnas[1]) > 0:
                plt.scatter(merged_df[args.x_col], merged_df[args.y_col],
                    edgecolor=sns.desaturate('black', 0.75), facecolors='none',
                    s=0.5)
                plt.scatter(snrnas[0], snrnas[1],
                            facecolor=facecolor, edgecolor=edgecolor,
                            label=sequence, marker=marker,
                            s=30, linewidth='1.5')

        if args.legend:
            l = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
            ncol=3, borderaxespad=0., fancybox=True,
            shadow=True, scatterpoints=1, fontsize=13)
            for text in l.get_texts():
                row = U_SNRNA_DATAFRAME[U_SNRNA_DATAFRAME.apply(
                    lambda x: x['SEQUENCE'] == text._text, axis=1)]
                text.set_color(row['COLOR'][row.index[0]])
    if args.color_threeprime_startswith is not None:
        sub_df = sub_df[merged['THREEPRIME_SEQUENCE_x'].map(
            lambda x: x.startswith(args.color_threeprime_startswith))]
        plt.scatter(sub_df[args.x_col], sub_df[args.y_col],
                    facecolor='r', edgecolor='r',
                    label='G downstream of cellular fragment')
    if args.spearman_r:
        ax.text(0.02, 0.97,
                (r'$\rho_s = {:5.2f}$'.format(
                    spearmanr(np.array(merged_df[args.x_col]),
                np.array(merged_df[args.y_col]))[0])),
                ha='left', va='center', transform=ax.transAxes,
                fontsize=14)
    if args.n_stats:
        ax.text(0.02, 0.925,
                '$n_x = {:.2E}$'.format(n_x).replace('+', ''),
                ha='left', va='center', transform=ax.transAxes,
                fontsize=14)
        ax.text(0.02, 0.88,
                '$n_y = {:.2E}$'.format(n_y).replace('+', ''),
                ha='left', va='center', transform=ax.transAxes,
                fontsize=14)
        ax.text(0.35, 0.97,
                '$u_x = {:.2E}$'.format(u_x).replace('+', ''),
                ha='left', va='center', transform=ax.transAxes,
                fontsize=14)
        ax.text(0.35, 0.925,
                '$u_y = {:.2E}$'.format(u_y).replace('+', ''),
                ha='left', va='center', transform=ax.transAxes,
                fontsize=14)
        ax.text(0.35, 0.88, '$u_{{merged}} = {:.2E}$'.format(len(merged_df)).replace('+', ''),
                ha='left', va='center', transform=ax.transAxes,
                fontsize=14)
    try:
        bbox_extra_artists = (l,)
    except NameError:
        bbox_extra_artists = None
    plt.savefig(args.outfile, format=args.format, dpi=args.dpi,
                bbox_extra_artists=bbox_extra_artists, bbox_inches='tight')
    plt.close()
    
def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-a', '--file-a', type=argparse.FileType('r'))
    parser.add_argument('-b', '--file-b', type=argparse.FileType('r'))
    parser.add_argument('-o', '--outfile', type=str, required=True)
    parser.add_argument('--on', type=str, default='SEQUENCE')
    parser.add_argument('--quantile-normalize', action='store_true')
    parser.add_argument('--x-col', type=str, default='COUNT')
    parser.add_argument('--y-col', type=str, default='COUNT')
    parser.add_argument('--how', type=str, default='inner')
    parser.add_argument('--dpi', type=int, default=None)
    parser.add_argument('--bandwidth', default=None, type=float)
    parser.add_argument('--color-threeprime-startswith', type=str,
                        default=None)
    parser.add_argument('--greater-than-zero')
    parser.add_argument('--add-histograms', action='store_true')
    parser.add_argument('--only-starts-with', type=str, nargs='+', default=None)
    parser.add_argument('--cutoff-cpm', type=float, default=None)
    parser.add_argument('--cutoff-reads', type=int, default=None)
    parser.add_argument('--no-header', action='store_true')
    parser.add_argument('--xlabel', type=str, default='')
    parser.add_argument('--ylabel', type=str, default='')
    parser.add_argument('--title', type=str, default='')
    parser.add_argument('--spearman-r', action='store_true')
    parser.add_argument('--min-val', type=int, default=1)
    parser.add_argument('--max-val', type=int, default=10000000)
    parser.add_argument('-f', '--format', default='eps')
    parser.add_argument('--highlight', action='store_true')
    parser.add_argument('--snRNAs', type=str, nargs='+')
    parser.add_argument('--uppercase-columns', action='store_true')
    parser.add_argument('--gzipped', action='store_true', help='If specified, unzips files first.')
    parser.add_argument('--legend', action='store_true', required=False)
    parser.add_argument('--plot-relatives', action='store_true')
    parser.add_argument('--n-stats', action='store_true')
    parser.add_argument('--restrict-lengths', nargs=2, type=int)
    parser.add_argument('--include-par', default=None, type=str, help='GCAAAAGCAG or GCGAAAGCAG')
    parser.add_argument('--cpm', action='store_true', help='Counts per million')
    parser.add_argument('--generate-eps', action='store_true')
    args = parser.parse_args(argv)
    process_files(args)

if __name__ == '__main__':
    main(sys.argv[1:])
