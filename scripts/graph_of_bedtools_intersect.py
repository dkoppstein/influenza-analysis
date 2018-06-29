import argparse
import sys
import pandas as pd
from pandas import DataFrame, Series
import matplotlib.pyplot as plt
from mpltools import style
from collections import defaultdict
from cStringIO import StringIO
from subprocess import Popen, PIPE

def mask(df, f):
    return df[f(df)]

def get_reads(from_table, corresponding_to):
    '''
    Gets all the reads from table "from_table" that have the same read
    name as in table "corresponding_to". 
    Returns a pandas DataFrame.
    '''
    from_table_set = set(from_table['READ_NAME'])
    corresponding_to_set = set(corresponding_to['READ_NAME'])
    reads_to_consider = from_table_set.intersection(corresponding_to_set)
    return pd.DataFrame(row for idx, row in from_table.iterrows()
                        if row['READ_NAME'] in reads_to_consider)

def limit_table(table, key, gt, lt):
    return table[(table[key] >= gt) & (table[key] <= lt)]

def restrict_table_by_length(table, length):
    return table[table['ORIGINAL_SEQUENCE'].apply(len) >= length]

def calculate_num_times_map(df):
    '''Adds a column, num_times_map, to the existing DataFrame.'''
    num_times_map = df.groupby('READ_NAME').size()
    df = df.set_index('READ_NAME')
    df['NUM_TIMES_MAP'] = num_times_map
    df['NUM_TIMES_MAP'][df['TSS_CHROM'] == '*'] = 0
    df = df.reset_index()
    return df

def graph_original_and_shuffled(d, args, length=None, num_times_map=None):
    '''Takes a dictionary containing the original and shuffled tables
    as values, as well as the desired length and/or
    num_times_map cutoff, and graphs it.'''
    d_copy = d.copy()
    fig, ax1 = plt.subplots()
    ax1.set_xlabel('Position relative to annotated TSS')
    if args.use_counts:
        ax1.set_ylabel('Number of counts')
    else:
        ax1.set_ylabel('Number of unique sequences')
    ax2 = ax1.twinx()
    ax2.set_ylim([0, 15])
    ax2.set_ylabel('Signal to noise', color='c')
    for ax in [ax1, ax2]:
        ax.set_xlim([-25,25])

    for label, table in d_copy.items():
        # restrict mapping to [-25, 25]
        table = limit_table(table, 'OFFSET_FROM_START', -25, 25)

        if num_times_map is not None:
            # first calculate num_times_map
            table = calculate_num_times_map(table)
            table = table[table['NUM_TIMES_MAP'] <= num_times_map]

        # restrict table by length
        if length is not None:
            table = restrict_table_by_length(table, length)
            
        # normalize each read by number of times it maps
        
        table['FREQUENCY'] = 1.
        if args.normalize_by_num_maps:
            table['FREQUENCY'] = table['FREQUENCY'] / \
              table['NUM_TIMES_MAP'].astype(float)

        if args.use_counts:
            table['FREQUENCY'] = table['FREQUENCY'] * table['COUNT']
            
        # sum the weighted reads and plot
        groupby_sum = table['FREQUENCY'].groupby(
            table['OFFSET_FROM_START']).sum()

        ax1.plot(groupby_sum.index, groupby_sum, label=label)

        # reassign the dictionary
        d_copy[label] = groupby_sum.reindex(range(-25, 26))

    signal_to_noise = d_copy['Original'] / d_copy['Shuffled']

    ax2.plot(d_copy['Original'].index, signal_to_noise,
        'c', label='Signal to noise')
    ax1.legend()
    plt.title(args.title)
    outfile = args.outfile[0]
    if length is not None:
        outfile += '_length_%d' % length
    if num_times_map is not None:
        outfile += '_numtimesmap_%d' % num_times_map
    outfile += '.%s' % args.format
    plt.savefig(outfile, format=args.format)

def process_file(args):
    style.use('ggplot')
    args.unshuffled
    if args.gzipped:
        args.unshuffled, args.shuffled = (StringIO(Popen(['zcat', f.name],
                                          stdout=PIPE).communicate()[0]) \
                                          for f in [args.unshuffled,
                                                    args.shuffled])
    original, shuffled = (pd.read_table(fn, index_col=False) for
                          fn in (args.unshuffled, args.shuffled))
    d = {'Original': original, 'Shuffled': shuffled}

    if args.lengths is not None:
        for length in range(*args.lengths):
            graph_original_and_shuffled(d, args, length=length)

    if args.num_times_map is not None:
        for num in range(*args.num_times_map):
            graph_original_and_shuffled(d, args, num_times_map=num)

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=\
                                     argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-u', '--unshuffled', type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-s', '--shuffled', type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-o', '--outfile', nargs=1, type=str,
                        required=True)
    parser.add_argument('--lengths', nargs=2, type=int)
    parser.add_argument('--num-times-map', nargs=2, type=int)
    parser.add_argument('-f', '--format', type=str,
                        default='png', nargs='?')
    parser.add_argument('-t', '--title', type=str, default='', nargs='?')
    parser.add_argument('--use-counts', default=False, action='store_true')
    parser.add_argument('--normalize-by-num-maps', default=False,
                        action='store_true')
    parser.add_argument('-m', '--max', type=int, default=20, nargs='?')
    parser.add_argument('--gzipped', action='store_true')
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
