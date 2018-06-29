#!/usr/bin/env python
# @author: David Koppstein

'''
Finds all the transcription start sites in a given gtf file. 

The way we determine starts is derived from:

http://www-huber.embl.de/users/anders/HTSeq/doc/tss.html

Briefly, we group by the "transcript id" feature, and then find the 
exon 

If -n is specified, adds n nucleotides to the end coordinates, and
subtracts n from the start coordinates. 
'''

import sys
import argparse
import os

# currently using HTSeq 0.5.4p3
import HTSeq

# types of gtf files to accept
# Refseq also refers to UCSC
# GTF_TYPES = ('Refseq', 'Ensembl')

def test():
    'Implement testing for this script at some point. Use nosetests.'

def process_refseq(args):
    '''
    Refseq files don't have exon_number attributes. 
    Therefore, we need to figure things out for ourselves.
    '''
    gtffile = HTSeq.GFF_Reader(args.infile[0])
    txs = {}
    for feature in gtffile:
        tid = feature.attr['transcript_id']
        if feature.type == 'exon':
            if tid not in txs:
                txs[tid] = [feature]
            else:
                txs[tid].append(feature)

    tsspos = {}
    for tid, lst in txs.iteritems():
        bestft = find_start(lst)
        # if it's on the negative strand, switch it (?)
        if bestft.iv.start_d_as_pos not in tsspos:
            start_start(bestft, args.num)
            tsspos[bestft.iv.start_d_as_pos] = bestft
    for ft in tsspos.itervalues():
        print ft.get_gff_line().rstrip()

def start_start(ft, num):
    '''
    Change the end position to the start position. Add num to the end, 
    subtract num from the beginning. If the beginning is less than 0, make it 0.
    '''
    if ft.iv.start == ft.iv.start_d:
        ft.iv.end = ft.iv.start+1+num
        ft.iv.start -= num
    else:
        ft.iv.start = ft.iv.end-1-num
        ft.iv.end += num
    if ft.iv.start < 0:
        ft.iv.start = 0

def find_start(lst):
    '''
    Takes a list of genomic features and returns the feature with the
    start coordinate, depending on strandedness.
    '''
    # where x is the old coordinate, y is the new. If comp_d[strand](x,y)
    # is true, then replace x with y. 
    comp_d = {'+': lambda x,y: x > y, '-': lambda x,y: x < y}

    bestft = None
    for ft in lst:
        if bestft is None:
            bestft = ft
        else:
            assert bestft.iv.strand == ft.iv.strand
        if comp_d[ft.iv.strand](bestft.iv.start_d, ft.iv.start_d):
            bestft = ft
    return bestft

# def process_ensembl(gtffile, num):
#     '''
#     Takes a GFF_Reader object, then prints output. Will be unsorted, 
#     so immediately pass to sort on command line. 
#     This is the way that is documented here:
#     http://www-huber.embl.de/users/anders/HTSeq/doc/tss.html
#     '''
#     tsspos = {}
#     for feature in gtffile:
#         if feature.type == 'exon' and feature.attr['exon_number'] == '1':
#                 if feature.iv.start_d_as_pos not in tsspos:
#                     tsspos[feature.iv.start_d_as_pos] = feature
#     for ft in tsspos.itervalues():
#         start_start(ft, num)
#         print ft.get_gff_line().rstrip()

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('infile', nargs=1, help='An ensembl, ' + \
                        'refseq, or ucsc gtf file')
    # parser.add_argument('-t', '--type', nargs=1, help='Specify one of ' + \
    #                     '"Refseq" (also includes UCSC) or "Ensembl"')
    parser.add_argument('-u', '--upstream', nargs='?', type=int, default=0, 
                        help='Number of nucleotides to add upstream of the TSS')
    parser.add_argument('-d', '--downstream', nargs='?', type=int, default=0,
                        help='Number of nucleotides to add downstream of the TSS')
    parser.add_argument('-n', '--num', nargs='?', type=int, default=0, 
                        help='Number of nucleotides to add to the ' + \
                        'coordinates. Deprecated as of Feb 2014.')
    # filetype = parser.parse_args().type[0]
    num = parser.parse_args().num

    # if filetype not in GTF_TYPES:
    #     raise RuntimeError("Filetype %s not supported, choose one of " + \
    #                        "%s" % (filetype, str(GTF_TYPES)))
    
    
    # process_func = dict(zip(GTF_TYPES, (process_refseq, process_ensembl)))
    # process_func[filetype](gtffile, num)
    args = parser.parse_args(argv)
    process_refseq(args)

if __name__ == "__main__":
    main(sys.argv[1:])
