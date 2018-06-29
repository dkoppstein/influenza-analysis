'''
Functions to turn sequences into PSSMs into stacked bar plots.

Also contains a command-line interface.
'''

__author__ = 'David Koppstein and Guo-Liang Chew'

import numpy as np
import pandas as pd
import sys
from flulib import fastq_to_dataframe
from pandas import DataFrame, Series
from scipy import log2
import matplotlib.pyplot as plt
import argparse
import matplotlib.patches as mpatches
import seaborn as sns
from cStringIO import StringIO
from subprocess import Popen, PIPE

FOURPEAKS_COLORS = {
    'A' : '#67f577',
    'C' : '#1f40ff',
    'G' : '#000000',
    'T' : '#f72e14',
    'N' : 'c'
    }

WEBLOGO_COLORS = {
    'A': '#00d700', # green
    'C': '#0226cc', # blue
    'G': '#ffb700', # yellow
    'T': '#df1f00' # red
    }

def normalize_lengths(iterable, shift_right=False):
    '''
    Given an iterable of sequences as strings, returns a pandas Series of sequences
    with Ns added the end to make them all the same length.
    
    If shift_right, adds Ns to the beginning instead of the end.
    '''
    series = Series(iterable)
    max_length = series.apply(len).max()
    def normalize_seq(seq):
        ns = 'N'*(max_length - len(seq))
        if shift_right:
            return ns + seq
        else:
            return seq + ns
    return series.apply(normalize_seq)
        
def explode_series(series, shift_right=False):
    '''
    Turns a pandas Series of sequences into a Dataframe,
    one nucleotide per position. 

    Conceptually, given two sequences AGTC and ACC, turns into: 
    
position   0 1 2 3
seq_name
seq1       A G T C
seq2       A C C N

    If shift_right is stipulated, returns: 

position   0 1 2 3
seq_name
seq1       A G T C
seq2       N A C C

    i.e. shortened sequences are normalized with respect to the three-prime
    end.
    '''
    normalized_series = normalize_lengths(series, shift_right=shift_right)
    return DataFrame(normalized_series.apply(
        lambda x: Series([letter.upper() for letter in x])))

def normalize_df(df):
    '''Normalizes a DataFrame by the total of each column.'''
    return df.div(df.sum(axis=0), axis=1)

def nucleotide_frequencies(df, count_series=None, normalize=False,
                           ignore_ns=False, weight_column='COUNT'):
    '''
    Gets the nucleotide frequencies, depicted as such:
    
position     0 1 2 3
letter
A            2 0 0 0
C            0 1 1 1
G            0 1 0 0
T            0 0 1 0
N            0 0 0 1

of an exploded series. If given a count_series, weights the ith sequence in
df by the ith value in count_series; otherwise, weights all sequences equally.

If normalize, will return the fraction for nucleotide instead:

position     0  1   2   3
letter
A            1. 0   0   0
C            0  0.5 0.5 0.5
G            0  0.5 0   0
T            0  0   0.5 0
N            0  0   0   0.5

    '''
    if count_series is not None:
        assert len(count_series) == len(df)
        df[weight_column] = count_series

    # get nucleotide frequencies per column
    def column_func(x): 
        if x.name == weight_column: # we're looking at the count column itself
            return None
        if count_series is None: # don't weight by counts
            return x.value_counts()
        else: # do weight by counts
            return df[weight_column].groupby(df[x.name]).sum()
        
    df = df.apply(column_func, axis=0).fillna(0)
    
    # remove the concatenated count series
    if weight_column in df.columns: 
        df.drop(weight_column, axis=1, inplace=True)

    if ignore_ns and 'N' in df.index:
        df.drop('N', axis=0, inplace=True)
        
    if normalize:
        df = normalize_df(df)
    return df

def plot_nucleotide_dist(pssm, outfile, format='png', title='',
                         max_y=0.7, weblogo_colors=True):
    sns.set_style("darkgrid", {"axes.facecolor": "#E5E5E5", 'font.family': 'Arial'})
    sns.set_context('paper')
    pssm = pssm.T
    if weblogo_colors:
        colors = WEBLOGO_COLORS
    else:
        colors = FOURPEAKS_COLORS
    fig, ax = plt.subplots()
    pssm.plot(color=[colors[i] for i in pssm.columns])
    plt.xlabel('Position')
    plt.ylabel('Frequency')
    plt.legend(loc='best', fontsize=16)
    plt.ylim(0.0, max_y)
    plt.title(title)
    plt.savefig(outfile, format=format)

def plot_dist_from_df(df, outfile, restrict_lengths=None,
                      format='png', title='', max_y=0.7, weblogo_colors=True,
                      shift_right=False, ignore_ns=False, ncols_from_right=None,
                      sequence_column='SEQUENCE', weight_column='COUNT'):
    sns.set_style("darkgrid", {"axes.facecolor": "#E5E5E5", 'font.family': 'Arial'})
    sns.set_context('paper')
    df.columns = [h.upper() for h in df.columns]
    if restrict_lengths:
        lower, upper = tuple(restrict_lengths)
        criterion = df[sequence_column].map(
            lambda x: len(x) >= lower and len(x) <= upper)
        df = df[criterion]
    series = explode_series(df[sequence_column], shift_right=shift_right)
    pssm = nucleotide_frequencies(series, count_series=df[weight_column],
                                normalize=True, ignore_ns=ignore_ns,
                                weight_column=weight_column)
    if ncols_from_right:
        pssm = pssm.ix[:, len(pssm.columns)-ncols_from_right:]
    plot_nucleotide_dist(pssm, outfile, format=format, title=title,
                         max_y=max_y, weblogo_colors=weblogo_colors)

def calc_log_likelihood(pssm, background=None, correct_lengths=False):
    if background is not None:
        if correct_lengths:
            background = background[range(background.columns[0],
                                          len(pssm.columns))]
            background.columns = pssm.columns
        assert type(background) is DataFrame
        if type(pssm.columns == background.columns) is bool:
            background = background[range(background.columns[0],
                                          len(pssm.columns))]
            background.columns = pssm.columns
        pssm /= background
    return pssm.apply(log2).replace(-np.inf, 0.)

def calculate_ic(pssm, background=None):
    '''Given a normalized PSSM, calculates the information content, given by

    IC(w) = log2(J) + \sum_{j=1}^{J} [p_{wj} * log2(p_{wj})]

    If given a background, computes

    IC(w) = log2(J) + \sum_{j=1}^{J} [p_{wj} * log2(p_{wj} / b_{wj})]
    '''
    ans = (pssm * calc_log_likelihood(pssm, background)).sum()
    if background is None:
        ans += log2(len(pssm))
    return ans

def create_background(df, infer=False):
    '''Creates a background distribution, with the same columns and
    rows as the original DataFrame.

    :param infer: If False, creates a uniform distribution
    (i.e. 0.25 per nucleotide). If True, creates a distribution with nucleotide
    frequencies equal to that in the given DataFrame. 
    :type infer: bool. 
    '''
    if infer:
        nuc_mean = df.mean(axis=1)
        background = DataFrame(nuc_mean for i in xrange(len(df.columns))).T
    else:
        nrows = len(df)
        background = DataFrame([[1. / nrows] * nrows] * len(df.columns)).T
    background.columns = df.columns
    background.index = df.index
    return background

def transform_pssm(pssm, background=None, normalize=False,
                   log_likelihood=False, ic_scale=True,
                   flat_background=False):
    '''
    Takes a PSSM and applies various transformations to it to get information content,
    log-likelihood, or a combination thereof. 

    Default values will give the frequency of the nucleotide in the PSSM scaled by
    the information content, but can be modified. 

    :param pssm: Position-weight matrix.
    :type pssm: DataFrame.
    :param background: Can be one of three things: None (uniform distribution;
    default), string 'infer' (calculates background distribution using
    nucleotide frequencies in PSSM), or a DataFrame. 
    :type background: NoneType, str, or DataFrame. 
    :param normalize: If specified, normalizes all columns to sum to one in the 
    input pssm and background. 
    :type normalize: bool. 
    :param log_likelihood: Instead of using the original values in the pssm to
    scale the nucleotide values, uses log2(p_{wj}/b_{wj}) where p is the
    probability of seeing nucleotide j at position w, and b is the corresponding
    probability of the background distribution.
    :type log_likelihood: bool.
    :param ic_scale: Scales each position by the information content, given by 

    IC(w) = log2(J) + \sum_{j=1}^{J} [p_{wj} * log2(p_{wj})]

    where w is the position, J is the number of nucleotides,
    and p_{wj} is the probability of seeing nucleotide j at position w.
    :type ic_scale: bool. 
    '''
    # create a background distribution where each nucleotide frequency
    # is the average across all positions
    if type(background) is str and background == 'infer':
        background = create_background(pssm, infer=True)
    if type(background) is DataFrame and flat_background:
        nuc_mean = background.mean(axis=1)
        background = DataFrame(nuc_mean for i in xrange(len(background.columns))).T
    if normalize: # normalize pssm and background
        try: 
            background, pssm = tuple(normalize_df(df) for df in (background, pssm))
        except NameError: # background is None
            pssm = normalize_df(pssm)

    if ic_scale: # calculate the information content before transforming the PSSM
        ic = calculate_ic(pssm, background)
        
    if log_likelihood: # log2(p/q)
        if background is None: # create a uniform distribution if need be
            background = create_background(pssm)
        pssm = calc_log_likelihood(pssm, background, correct_lengths=True)

    try: # scale the pssm by the precomputed information content
        pssm *= ic
    except NameError: # ic_scale was False
        pass
        
    return pssm

def round_to(x, prec=2, base=.05):
    ''''Courtesy of
    http://stackoverflow.com/questions/2272149/round-to-5-or-other-number-in-python
    '''
    return round(base * round(float(x) / base), prec)

def plot_stacked_bar(pssm, outfile, 
                     ylabel='Information Content (bits)', format='png',
                     max_x=None, use_bottom=True, force_yticks=False,
                     max_y=None, label_xticks=None,
                     offset_bars=False):
    '''
    :param pssm: Plots whatever is in pssm. Can modify with the various
    options in the above function, entropy. 
    '''
    sns.set_style("darkgrid", {"axes.facecolor": "#E5E5E5", 'font.family': 'Arial'})
    sns.set_context('paper')
    fig = plt.figure()
    ax = plt.gca()

    if max_x is not None:
        pssm = pssm.ix[:, :list(pssm.columns).index(max_x)]

    npos = len(pssm.columns) # number of positions to plot
    width = 0.2
    index = np.arange(0, width * npos, width)
    
    # no more than 5 colors (?)
    colors = [WEBLOGO_COLORS[x] for x in ('A', 'C', 'G', 'T')]
    colors.append('c') # for N
    colors = colors[:len(pssm.index)] 

    # use bottom of previous positive values for positive values,
    # and bottom of previous negative values for negative values
    max_val= 0.
    use_neg = False        
    for idx, i in enumerate(pssm.columns):
        # sort by decreasing absolute value (for visualization)
        temp_pssm = pssm.copy()
        temp_pssm['TEMP_COLUMN'] = pssm[i].map(abs)
        if not offset_bars: # if we offset, don't sort
            temp_pssm = temp_pssm.sort(columns=['TEMP_COLUMN'],
                                       ascending=False)
        pos_totals, neg_totals = 0., 0.
        pos_bottom, neg_bottom = (np.array([0.]*len(pssm.columns)) for i in range(2))
        for j, nuc in enumerate(temp_pssm.index):
            color = WEBLOGO_COLORS[nuc]
            val = temp_pssm[i].ix[nuc]
            if val < 0:
                use_neg = True
            abs_val = temp_pssm['TEMP_COLUMN'].ix[nuc]
            if abs_val > max_val:
                max_val = abs_val
            arr = np.array([0.]*len(pssm.columns))
            arr[idx] = val
            if val >= 0:
                bottom = pos_bottom
            else:
                bottom = neg_bottom
            if offset_bars:
                index_to_plot = np.array([(x*5.+j*width) for x in index])
            else:
                index_to_plot = index
            plt.bar(index_to_plot, arr, width, color=color, bottom=bottom,
                    label=nuc, edgecolor='none')
            if use_bottom:
                if val >= 0:
                    pos_bottom += arr
                else:
                    neg_bottom += arr

    # choose the max of the pos_totals and neg_totals
    if use_neg:
        min_val = -max_val
    else:
        min_val = 0.
    plt.ylabel(ylabel, fontsize=16)
    plt.xlabel('Position', fontsize=16)
    
    labels = ['A', 'C', 'G', 'T']
    patches = [mpatches.Patch(color=WEBLOGO_COLORS[k], label=k) for k \
               in labels]
    plt.legend(patches, labels, loc='best', fontsize=16)
    if max_y:
        max_val, min_val = max_y, -max_y
    plt.ylim(min_val, max_val)
    if label_xticks:
        xtick_labels = range(label_xticks[0], label_xticks[1])
    else:
        xtick_labels = [str(x) for x in pssm.columns]
    if offset_bars:
        plt.xticks(np.array([(x*5. + 9.*5.*width/10.) for x in index]),
                   xtick_labels)
    else:
        plt.xticks(index + width / 2., xtick_labels)
    # sensible number of ticks: about 20 per plot
    ytick_width = round_to((max_val - min_val) / 20., base=0.05)
    if ytick_width == 0:
        ytick_width += 0.05
    
    # create a range of yticks s.t. one of the ticks intersects 0.
    if min_val < 0:
        yticks = np.arange(0, min_val, -ytick_width)[::-1]
    else:
        yticks = np.arange(0)
    yticks = np.append(yticks, np.arange(0, max_val, ytick_width))
    if force_yticks:
        plt.yticks(np.arange(0, 0.05, 0.02))
    else:
        plt.yticks(yticks)

    for axis in (ax.xaxis, ax.yaxis):
        for tick in axis.get_major_ticks():
            tick.label.set_fontsize(16)
    plt.tight_layout()
    plt.savefig(outfile, format=format)
    plt.close()

def parse_pssm(infile):
    pssm = pd.read_table(infile, index_col=0) # first column is index
    pssm.columns = [int(x) for x in pssm.columns] # columns are ints
    return pssm

def parse_table(infile, restrict_lengths=None, ncols_from_right=None,
                      ignore_ns=True, no_plot_ns=False, shift_right=False,
                      restrict_end_to=None, sequence_column='SEQUENCE',
                      weight_column='WEIGHT', gzipped=False,
                      infer_columns=False):
    if type(infile) is not pd.DataFrame:
        if gzipped:
            infile = StringIO(Popen(['zcat', infile.name],
                            stdout=PIPE).communicate()[0])
        table = pd.read_table(infile, index_col=None)
    else:
        table = infile
    if type(table.columns[0]) is str:
        if infer_columns:
            if type(table[table.columns[0]][0]) is str:
                table.columns = [sequence_column, weight_column]
            else:
                table.columns = [weight_column, sequence_column]
        else:
            table.columns = [c.upper() for c in table.columns]
    elif len(table.columns) == 2:
        table.columns = [sequence_column, weight_column]
    else:
        raise 'Columns of input table must either be str or length 2.x'
    # restrict to not null
    table = table[table[sequence_column].notnull()]
    if restrict_lengths:
        lower, upper = tuple(restrict_lengths)
        criterion = table[sequence_column].map(lambda x: len(x) >= lower and len(x) <= upper)
        table = table[criterion]
    if restrict_end_to is not None:
        table = table[table['3PTRIMMED'] == restrict_end_to]
    series = explode_series(table[sequence_column], shift_right=shift_right)
    pssm = nucleotide_frequencies(series, count_series=table[weight_column],
                                normalize=True, ignore_ns=ignore_ns)
    if no_plot_ns and 'N' in pssm.index:
        pssm = pssm.drop('N')
    return pssm

def process_files(args):
    shift_right = args.ncols_from_right is not None
    if args.table:
        table_parser = lambda x: parse_table(x, args.restrict_lengths,
                                            args.ncols_from_right, args.ignore_ns,
                                            no_plot_ns=args.no_plot_ns,
                                            restrict_end_to=args.restrict_end_to,
                                            sequence_column=args.sequence_column,
                                            shift_right=shift_right,
                                            weight_column=args.weight_column)
    elif args.combine_left_and_right:
        def combiner(infile):
            left_pssm = parse_table(infile, args.restrict_lengths,
                                    args.ncols_from_right, args.ignore_ns,
                                    no_plot_ns=args.no_plot_ns,
                                    restrict_end_to=args.restrict_end_to,
                                    sequence_column='ORIGINAL_SEQUENCE',
                                    shift_right=shift_right,
                                    weight_column=args.weight_column)
            right_pssm = parse_table(infile, args.restrict_lengths,
                                     args.ncols_from_right, args.ignore_ns,
                                     no_plot_ns=args.no_plot_ns,
                                     restrict_end_to=args.restrict_end_to,
                                     sequence_column='THREEPRIME_OF_CLEAVAGE',
                                     shift_right=False,
                                     weight_column=args.weight_column)
            
            return left_pssm
        table_parser = combiner
    else:
        table_parser = parse_pssm
    pssm = table_parser(args.infile)
    if args.ncols_from_right:
        pssm = pssm.ix[:, (len(pssm.columns)-args.ncols_from_right):]
    if args.ncols_from_left:
        pssm = pssm.ix[:, :args.ncols_from_left]
    if args.negative_offsets:
        pssm.columns = range(-len(pssm.columns), 0)
    if args.plot_nucleotide_dist:
        plot_nucleotide_dist(pssm, args.outfile, format=args.format,
                             max_y = args.max_y)
        return
    if args.background_file is not None:
        if args.simple_background_parse:
            background = parse_table(args.background_file,
                                     gzipped=args.gzipped_background,
                                     weight_column=args.weight_column,
                                     infer_columns=True)
        else:
            background = table_parser(args.background_file)
    elif args.infer_background:
        background = 'infer'
    else:
        background = None
    pssm = transform_pssm(pssm, background=background, normalize=args.normalize,
                          log_likelihood=args.log_likelihood,
                          ic_scale=(not args.no_ic_scale),
                          flat_background=args.flat_background) # double negatives! 
    plot_stacked_bar(pssm, args.outfile, ylabel=args.ylabel, max_x=args.max_x,
                     use_bottom=not args.log_likelihood, format=args.format,
                     force_yticks=args.force_yticks, max_y=args.max_y,
                     label_xticks=args.label_xticks, offset_bars=args.offset_bars)

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin, help='PSSM.')
    parser.add_argument('--table', action='store_true')
    parser.add_argument('--plot-nucleotide-dist', action='store_true')
    parser.add_argument('--ignore-ns', action='store_true')
    parser.add_argument('-b', '--background-file', type=argparse.FileType('r'),
                        default=None, help=('Can also specify a PSSM of a '
                        'background distribution.'))
    parser.add_argument('--infer-background', action='store_true',
                        help=('Will infer background from nucleotide frequencies'
                        'in the given infile.'))
    parser.add_argument('--ncols-from-right', type=int, default=None)
    parser.add_argument('--sequence-column', type=str, default='SEQUENCE')
    parser.add_argument('--weight-column', type=str, default='COUNT')
    parser.add_argument('--ncols-from-left', type=int)
    parser.add_argument('--restrict-end-to', type=str)
    parser.add_argument('--combine-left-and-right', action='store_true')
    parser.add_argument('--normalize', action='store_true')
    parser.add_argument('--restrict-lengths', type=int, nargs=2)
    parser.add_argument('--log-likelihood', action='store_true')
    parser.add_argument('--ylabel', type=str, default='Information Content (bits)')
    parser.add_argument('--max-x', type=int, nargs='?')
    parser.add_argument('--force-yticks', action='store_true')
    parser.add_argument('--max-y', type=float, nargs='?')
    parser.add_argument('--no-plot-ns', action='store_true',
                        help=('Unlike --ignore-ns, --no-plot-ns removes the Ns *after*'
                        'normalization, not before.'))
    parser.add_argument('--no-ic-scale', action='store_true')
    parser.add_argument('-o', '--outfile', type=str)
    parser.add_argument('--negative-offsets', action='store_true',
                        help='If true, uses negative offsets for the positions')
    parser.add_argument('--label-xticks', type=int, nargs=2, default=None)
    parser.add_argument('-f', '--format', type=str)
    parser.add_argument('--berrylogo', action='store_true')
    parser.add_argument('--offset-bars', action='store_true')
    parser.add_argument('--flat-background', action='store_true',
                        help='If true, average the background PSSM first.')
    parser.add_argument('--simple-background-parse', action='store_true',
                        help='Otherwise, defaults to using the '
                        'settings specified')
    parser.add_argument('--gzipped-background', action='store_true')
    parser.add_argument('--background-weight-column', type=str, default='COUNT')
    args = parser.parse_args(argv)
    process_files(args)
    
if __name__ == '__main__':
    main(sys.argv[1:]) 
