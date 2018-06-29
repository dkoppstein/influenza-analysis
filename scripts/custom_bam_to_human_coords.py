#!/usr/bin/env python
# custom_bam_to_human_coords.py
__author__ = 'David Koppstein'

'''
Converts the coordinates of a given BAM file from custom bowtie coordinates
to human coordinates. 
'''

import argparse
import sys
import HTSeq
import itertools
import re
from subprocess import Popen, PIPE
import sys

# SAM code for reverse complemented sequence
REVERSE_COMPLEMENT_FLAG = 16

def append_dist(aln):
    mapped_iv = parse_bowtie_interval(aln.iv.chrom)
    dist = aln.iv.start - int(mapped_iv.length / 2)
    aln.read.name += ':distance=' + str(dist)
    return

def reheader(genome_location):
    '''
    Creates a sam to bam compressing subprocess that has the human genome
    header information preloaded into it.
    '''
    cmd = ['bowtie', '-S', genome_location, '-m', '1', '-c', 'G'*1000]
    samtools = Popen(['samtools', 'view', '-bS', '-'],
                     stdin=PIPE, stdout=sys.stdout)
    
    # get the headers and pipe them into samtools view
    for header in Popen(cmd, stdout=PIPE).communicate()[0].splitlines():
        if header.startswith('@'):
            samtools.stdin.write(header + '\n')
        else:
            break
        
    return samtools

def parse_bowtie_interval(s):
    '''
    Returns a new HTSeq.GenomicInterval object based on the bowtie-formatted
    string given. 
    
    For example, 'chr7:5570181-5570282(-)' -> chr7:[5570181,5570282)/-
    '''
    pattern = '^(chr[1-9]|chr1[0-9]|chr2[0-2]|chrX|chrY|chrM):(\d+)-(\d+)' + \
    '\(([-+])\)$'
    match = re.search(pattern, s)
    if match is not None:
        chrom, start, end, strand = tuple(match.groups())
        return HTSeq.GenomicInterval(chrom, int(start), int(end), strand)
    else:
        return None
    
def process_file(infile, genome_location, distances):
    '''
    Converts the mapped BAM output of bowtie mapping to a custom bowtie
    index. 
    
    Calls bowtie on the fly to generate headers.

    Genome_location must be the path to a set of .ebwt bowtie indices whose
    genomic locations will be used for the new .bam header. 
    '''

    # make up a fake mapping to get the headers (can also use Samtools reheader?)
    samtools = reheader(genome_location)
    bam_reader = HTSeq.BAM_Reader(infile)
    for aln in bam_reader:
        # if not aligned, we can just write it as is
        if not aln.aligned:
            samtools.stdin.write(aln.get_sam_line() + '\n')
            continue
        if distances:
            append_dist(aln)
        mapped_iv = parse_bowtie_interval(aln.iv.chrom)
        if mapped_iv is None:
            continue
        if mapped_iv.strand == '+':
            s, e = mapped_iv.start + aln.iv.start, mapped_iv.end + aln.iv.end
        else:
            s, e = mapped_iv.end - aln.iv.start - aln.iv.length, mapped_iv.end - aln.iv.start + 1
        new_iv = HTSeq.GenomicInterval(mapped_iv.chrom, s, e, mapped_iv.strand)
        aln.iv = new_iv
        # get_sam_line doesn't set the flag according to strandedness
        # so we need to do things ourselves
        if new_iv.strand == '-':
            aln.flag += REVERSE_COMPLEMENT_FLAG
        samtools.stdin.write(aln.get_sam_line() + '\n')
        
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('infile', nargs=1, help='Input file.')
    parser.add_argument('-g', '--genome', nargs=1, help=('Human genome basename.'
                        ' Must be a path to a bowtie ebwt basename.'))
    parser.add_argument('-d', '--distances', action='store_true', 
                        help='Add distance from TSS for each alignment.')
    args = parser.parse_args()
    process_file(args.infile[0], args.genome[0], args.distances)

if __name__ == '__main__':
    main()
