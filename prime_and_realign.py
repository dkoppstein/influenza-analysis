from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import itertools
import sys
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq
from matplotlib.font_manager import FontProperties
from scipy.stats import rankdata
import matplotlib.colors as mcolors
from matplotlib import cm

# colors

def make_colormap(first_color, second_color):
    """Return a LinearSegmentedColormap
    seq: a sequence of floats and RGB-tuples. The floats should be increasing
    and in the interval (0,1).

    From http://stackoverflow.com/questions/16834861/create-own-colormap-using-matplotlib-and-plot-color-scale
    """
    seq = [first_color, second_color]
    seq = [(None,) * 3, 0.0] + list(seq) + [1.0, (None,) * 3]
    cdict = {'red': [], 'green': [], 'blue': []}
    for i, item in enumerate(seq):
        if isinstance(item, float):
            r1, g1, b1 = seq[i - 1]
            r2, g2, b2 = seq[i + 1]
            cdict['red'].append([item, r1, r2])
            cdict['green'].append([item, g1, g2])
            cdict['blue'].append([item, b1, b2])
    return mcolors.LinearSegmentedColormap('CustomMap', cdict)

# Functions for dataframes
def sort_and_rank(df, nrows=4):
    '''
    Add the rank per group, and restricts the groupby. 
    '''
    df = df.sort('COUNT', ascending=False)
    df['RANK'] = [x+1 for x in range(len(df))]
    return df[:nrows]

class NucNode:
    ''''A node in a tree of DNA sequences. Each node corresponds to a particular
    sequence, but contains only a single nucleotide.  

    For example: 

      A 
   /    \
  C      T
 / \    / \
A   T  G   T

The "A" on the bottom left corresponds to the sequence "ACA" 
(read from the top down)

This tree structure enables finding prime-and-realign events. 
'''

    def __init__(self):
        # By default, the children become a new NucNode
        self.children = defaultdict(NucNode)
        self.count = 0
        self.nuc = None
        self.parent = None
        self.depth = 0

    def update(self, sequence, count):
        sequence = sequence.upper()
        if sequence == '':
            self.count += count
        else:
            # Give yourself the first nucleotide
            nuc = sequence[:1] 
            # Recursively call this function on the children with the rest of the
            # nucleotides
            self.children[nuc].update(sequence[1:], count) 
            # Update the children
            self.children[nuc].nuc = nuc
            self.children[nuc].parent = self
            self.children[nuc].depth = self.depth + 1

    def __str__(self):
        return 'seq: ' + self.get_sequence() + '; count: ' + str(self.count) + \
          '; Children: ' + ''.join([str(child) for child in 
                                    self.children.itervalues()])

    def endswith(self, seq):
        return self.get_sequence().endswith(seq)
    
    def find_gcfs(self, cutoff=1, has_children=False, with_end=None,
                  gcf_not_end_with=None):
        '''
        Returns a list of pointers to NucNodes that are the "greatest common
        factors", i.e. they are the first ones to contain a count at that level.

        If cutoff is specified, requires at least that many counts.

        If has_children is True, requires the GCF to have children to be output.

        If with_end is not None, requires the GCF to have a child
        that has a particular end (for example, GCA).

        If at_depth is specified, then requires the GCF to have the child with
        a particular sequence at a certain depth (i.e. GCA at depth 3). 

        If gcf_not_end_with is specified, the GCF itself cannot end with a
        particular nucleotide (i.e. GCA). 
        '''
        if with_end is not None:
            at_depth = len(with_end)
        lst = []
        self._gcf_helper(lst, cutoff=cutoff, has_children=has_children),
        if with_end is not None:
            lst = [node for node in lst if \
                   node.has_child_with_end(with_end, at_depth=at_depth,
                                           cutoff=cutoff)]
        if gcf_not_end_with is not None:
            lst = [node for node in lst if \
                   not node.get_sequence().endswith(gcf_not_end_with)]
        return lst

    def has_child_with_end(self, seq, cutoff=1, at_depth=None):
        lst = [child for child in self.all_children() if \
                child.endswith(seq) and child.count >= cutoff]
        if at_depth is not None:
            lst = [child for child in lst if child.depth - self.depth == at_depth]
        return len(lst) > 0
    
    def _gcf_helper(self, lst, cutoff=1, has_children=False,
                    with_end=None, at_depth=None):
        'Helper function for find_gcfs. '
        if self.count >= cutoff:
            if not has_children or len(self.all_children()) > 0:
                lst.append(self)
        else: # recurse
            for n in self.children.itervalues():
                n._gcf_helper(lst, cutoff=cutoff, has_children=has_children)

    def get_children(self):
        return self.children.values()
                
    def get_sequence(self):
        node, seq = self, ''
        while node.nuc is not None:
            seq, node = node.nuc + seq, node.parent # walk down the tree
        return seq

    def all_children(self):
        'Return a list of *all* the children below this NucNode.'
        return self._create_lst(self._all_children_helper)

    def _create_lst(self, func):
        lst = []
        func(lst)
        return lst
    
    def _all_children_helper(self, lst):
        for n in self.children.itervalues():
            lst.append(n)
            n._all_children_helper(lst)

    def print_node(outh):
        pass
            
    def populate_node(self, target):
        '''
        Populate the target node with all the information at and below this node.
        '''
        target.count += self.count
        length = len(self.get_sequence())
        for child in self.all_children():
            target.update(child.get_sequence()[length:], child.count)

    def max_length(self):
        '''
        Return the length of the longest sequence below this NucNode.
        '''
        best = 0
        for child in self.all_children():
            length = len(child.get_sequence())
            if length > best: best = length
        return best

    def get_max_counts(self):
        df = self.get_dataframe()
        return max(df.max())
    
    def get_dataframe(self, nrows=4, normalize_to_self=True):
        '''Turns a nucnode into a pandas dataframe:
        columns: sequence, length, count, rank
        :param normalize_to_self: If specified, will normalize the count of each
        of the children of the node to the count of this node.
        :type normalize_to_self: bool.
        '''
        d = defaultdict(dict)
        for index, child in enumerate(self.all_children()):
            seq = child.get_sequence()
            d['SEQUENCE'][index] = seq
            d['LENGTH'][index] = len(seq)
            d['COUNT'][index] = child.count
        df = pd.DataFrame.from_dict(d)
        df = df[df['COUNT'] > 0] # discard things without counts
        if normalize_to_self:
            df['COUNT'] = df['COUNT'] / self.count

        # calculate rank of count within 
        return df.groupby('LENGTH').apply(
            lambda x: sort_and_rank(x, nrows=nrows))

    def has_children(self):
        return len(self.all_children()) > 0

    def write_seqrecord(self, outh, fmt='fasta'):
        # Just make up quality scores for now if FASTQ
        for i in range(self.count): 
            record = SeqRecord(Seq(self.get_sequence(), IUPACAmbiguousDNA()))
            if fmt.startswith('fastq'):
                record.letter_annotations = {'phred_quality': [10]*len(record)}
            outh.write(record.format(fmt))

    def heatmap_with_labels(self, outfile, text_color='black', edgecolors='w',
                            nrows=4, lengths=[2, 8], color_map='Pastel1',
                            format='eps'):
        '''Modified from here:
        http://stackoverflow.com/questions/21024066/annotate-heatmap-with-value-from-pandas-dataframe

        Used for taking the aggregated "stats" node and creating a heatmap
        of all the sequences below it. 
        '''
        df = self.get_dataframe(nrows=nrows) # sequence, length, count, rank
        df = df[(df['LENGTH'] >= lengths[0]) & (df['LENGTH'] <= lengths[1])]
        # reshape s.t. index: rank; columns: length.
        # can select subsets with df['count'] and df['sequence']
        df = df.pivot(index='RANK', columns='LENGTH')
        df = df.sort_index(ascending=False)
        width = len(df.columns)/7*10
        height = len(df.index)/7*10
    
        fig, ax = plt.subplots(figsize=(20,10)) # (figsize=(width,height))

        # normalize count to [0, 1]
        count_df = df['COUNT'].copy().fillna(0.).astype(float)

        # make color map
        c = mcolors.ColorConverter().to_rgb
        # cmap = make_colormap(c('#0000FF'), c('#FFFF00')) # blue-yellow
        cmap = plt.get_cmap(color_map)

        # put white lines between squares in heatmap
        heatmap = ax.pcolor(count_df, cmap=cmap)

        for x_idx, x in enumerate(df['SEQUENCE'].index):
            for y_idx, y in enumerate(df['SEQUENCE'].columns):
                plt.text(y_idx + 0.5, x_idx + 0.5,
                         '%s' % str(df['SEQUENCE'][y][x]), 
                          horizontalalignment='center', fontsize=20,
                          verticalalignment='center', color=text_color)

        ax.autoscale(tight=True)  # get rid of whitespace in margins of heatmap
        ax.set_aspect('equal')  # ensure heatmap cells are square
        ax.xaxis.set_ticks_position('top')  # put column labels at the top
        # turn off ticks
        ax.tick_params(bottom='off', top='off', left='off', right='off')

        ax.set_yticks(np.arange(len(count_df.index)) + 0.5)
        ax.set_yticklabels(count_df.index, size=20)
        ax.set_xticks(np.arange(len(count_df.columns)) + 0.5)
        ax.set_xticklabels(count_df.columns, size=24)
        plt.xlabel('Nucleotides Realigned', fontsize=36)
        ax.xaxis.set_label_position('top')
        plt.ylabel('Rank Within Group', fontsize=36)
        
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", "3%", pad="1%")
        cbar = plt.colorbar(heatmap, cax=cax)
        cbar.set_label(r'Fraction realigned', fontsize=36)
        for t in cbar.ax.get_yticklabels():
            t.set_fontsize(20)
        plt.savefig(outfile, format=format)


def list_of_nucs(nucs, length):
    '''
    Creates a list of all possible nucleotides at or less than a given length
    >>> list_of_nucs('ACGT', 2)
    ['A', 'C', 'G', 'T', 'AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 
    'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT']
    '''
    lst = []
    for i in range(1, length+1):
        lst.extend(itertools.product(nucs, repeat=i))
    return [''.join(x) for x in lst]


def dataframe_from_gcfs(node, args):
    'Gets sequence and count'
    d = defaultdict(dict)
    for i, gcf in enumerate(node.find_gcfs(cutoff=args.cutoff,
                            has_children=args.has_children,
                            with_end=args.with_end,
                            gcf_not_end_with=args.gcf_not_end_with)):
        d['SEQUENCE'][i] = gcf.get_sequence()
        d['COUNT'][i] = gcf.count
    df = pd.DataFrame.from_dict(d)
    df = df.sort_index(ascending=False, axis=1) # sequence, then count
    return df


def table_from_seqrecords(infile, fmt): 
    l = []
    for record in SeqIO.parse(infile, fmt):
        l.append(str(record.seq))
    table = pd.DataFrame(pd.Series(l).value_counts()).reset_index()
    table.columns = ['SEQUENCE', 'COUNT']
    return table


def process_file(args):
    # modified to use full_table.txt
    # base node to add sequences to
    base = NucNode()
    table = table_from_seqrecords(args.infile, args.input_format)

    # get the original sequence, minus any 5-prime-end trimmed Gs
    for sequence, count in zip(table['SEQUENCE'], table['COUNT']):
        base.update(sequence, count) # populate the tree
    total_reads = table['COUNT'].sum() # do we need this? 

    # make a new node to aggregate information from GCFs in base node
    stats_node = NucNode()
    for node in base.find_gcfs(cutoff=args.cutoff,
                               has_children=args.has_children,
                               with_end=args.with_end):
        node.populate_node(stats_node)

    stats_node.heatmap_with_labels(args.outfile,
                                   nrows=args.nrows,
                                   lengths=args.lengths,
                                   color_map=args.colormap,
                                   format=args.output_format)


def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)

    parser.add_argument('--use-counts', action='store_true')
    parser.add_argument('--cutoff', type=int, nargs='?', default=1,
                        help='A nuc node must have >= cutoff reads to be ' + \
                        'considered a GCF node.')
    parser.add_argument('-o', '--outfile', type=str,
                        required=True)
    parser.add_argument('--as-table', action='store_true')
    parser.add_argument('--output-gcfs', action='store_true',
                        help='Write a pandas table the GCFs')
    parser.add_argument('--as-lengths', action='store_true',
                        help='Output length distribution of GCFs.')
    parser.add_argument('--as-pssm', action='store_true')
    parser.add_argument('--output-format', default='png', type=str, nargs='?',
                        help='Format of output file.')
    parser.add_argument('--input-format', default='fastq', type=str, nargs='?',
                        help='Format of input file.')
    parser.add_argument('--has-children', action='store_true',
                        help='Requires each GCF to have at least one child.')
    parser.add_argument('--with-end', type=str, nargs='?', required=False)
    parser.add_argument('--colormap', type=str, default='coolwarm')
    parser.add_argument('--gcf-not-end-with', type=str, nargs='?', required=False)
    parser.add_argument('--last-nucs', type=int, nargs='?',
                        help='Last n nucleotides to output GCFs.')
    parser.add_argument('--nrows', type=int, default=5,
                        help=('Number of sequences per prime-and-realign '
                                'length'))
    parser.add_argument('--lengths', type=int, nargs=2, default=[1, 8],
                        help='Number of lengths to graph when plotting heatmap.')
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
