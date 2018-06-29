from cStringIO import StringIO
from subprocess import Popen, PIPE
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import Bio.SeqIO
import sys
from collections import defaultdict
from mpltools import style

def plot_length_dist(series, fmt='eps', axis=None,
                     count_series=None, xlims=[0,20],
                     title='', label='', color='k'):
    'Takes a series of sequences, and plots their lengths.'
    style.use('ggplot')
    series = series.apply(len)
    if count_series is None: 
        count_series = series.apply(lambda x: 1)
    length_df = pd.concat([series, count_series], axis=1)
    length_df.columns = ['Length', 'Count']
    final_series = length_df.groupby('Length').sum()
    if 0 not in final_series.index:
        final_series.ix[0] = 0
    final_series = final_series.sort()
    final_series.columns = [label]
    ax = final_series.plot(xticks=range(*xlims), color=color, ax=axis)
    plt.xlim(*xlims)
    if len(title) > 0:
        plt.title(title)
    plt.xlabel('Cellular Fragment Length', fontsize=16)
    # get the exponent and label it
    ax = plt.gca()
    ax.get_yaxis().get_major_formatter().set_scientific(True)
    offset = ax.get_yaxis().get_offset_text()
    plt.ylabel(r'Counts ( x $10^7$)', fontsize=16)
    return ax

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=\
                                     argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infiles-untrimmed', nargs='+',
                        type=argparse.FileType('r'), help='Input files.')
    parser.add_argument('--infiles-trimmed', nargs='+',
                        type=argparse.FileType('r'),
                        help='Input files after prime-and-realign trimming.')
    parser.add_argument('--gzipped', action='store_true')
    parser.add_argument('-f', '--format', nargs='?', type=str,
                        default='fastq-illumina', 
                        choices=['fastq-illumina', 'fastq-solexa', 'fastq',
                                 'table', 'fasta'],
                        help='Format of input file.')
    parser.add_argument('-o', '--outfile', nargs='?', 
                        type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('-u', '--output-format', type=str, default='eps',
                        nargs='?')
    parser.add_argument('-t', '--title', type=str, default='', 
                        nargs='?')
    args = parser.parse_args()
    ax = None
    for legend_name, color, file_type in \
        [('Before trimming', 'black', args.infiles_untrimmed),
         ('After trimming', 'blue', args.infiles_trimmed)]:
        series, count_s = pd.Series([]), pd.Series([])
        for f in file_type:
            if args.gzipped:
                f = StringIO(Popen(['zcat', f.name],
                            stdout=PIPE).communicate()[0])
            if args.format == 'table':
                df = pd.read_table(f, index_col=None)
                df.columns = [str(s).upper() for s in df.columns]
                if 'SEQUENCE' not in df.columns:
                    df.columns = ['SEQUENCE', 'COUNT']
                series = pd.concat([series, df['SEQUENCE']])
                count_s = pd.concat([count_s, df['COUNT']])
            else:
                d_dict = defaultdict(int)
                for record in Bio.SeqIO.parse(f, args.format):
                    d_dict[len(record)] += 1
                if not 0 in d_dict:
                    d_dict[0] = 0
                series = pd.concat([series, pd.Series(d_dict)])
                count_s = None
        ax = plot_length_dist(series, fmt=args.output_format,
                            count_series=count_s, title=args.title, label=legend_name, color=color, axis=ax)
    ax.legend(loc='best')
    plt.savefig(args.outfile, format=args.output_format)

if __name__ == '__main__':
    main()
