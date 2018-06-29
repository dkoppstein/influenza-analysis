#!/usr/bin/env python
# count_transcripts.py

__author__ = 'David Koppstein'

r'''
Input: a fastq or fasta file.

Uses the first "nucs" nucleotides of the read as a barcode to count
cDNA molecules. If there are identical reads with the same barcode, only the
highest quality one is kept. Outputs a collapsed file that is missing the
barcodes.

If --guanines is stipulated, will additionally trim guanines after the barcode
and use them for barcode counting as well.

>> main(['test/test_multiple_in.fastq'])
@HWI-ST978_0130:2:1101:1496:2431#0/6
GCAGCACGTAAATATTGGCGAAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACAGATCTCGTATGCCGTCTTCTGCTT
+
eecgehhhfdegfh]gghhhhhcfhhhhh`db_Z]baa[_FGTX_aaaabRYRW[^BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
@HWI-ST978_0130:2:1101:1435:2448#0/7
GCTTATCAGACTGATGTTGACAAAAAAAAAAAAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGAACGATCTCGTATGCCGTCTTTT
+
gegfghiiihiiiiiiiiiiiiiiiihihhiiedZ]V^`Z`aZ^\^a\aW_bbY^`^^]``ba_RTTY]EJGW`[^O]b_bcacacBB
'''

from collections import defaultdict
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import numpy as np
import re
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.Seq import Seq

def get_count(record_id):
    '''Return an int of the count in the record header. DEPRECATED; use flulib instead.
    Or whatever I'm going to call it.
    '''
    pat = r'count=(\d+)'
    return int(re.search(pat, record_id).groups()[0])

def mean_quality(record):
    'Return the mean Phred quality of the given SeqRecord.'
    return np.mean(record.letter_annotations['phred_quality'])

def process_file(args):
    if args.barcode_file is None and args.barcode_sequence is not None:
        parser.error('If specifying a sequence to output barcodes for, you ' + \
                     'need to specify a --barcode-file as well!')

    if args.make_table:
        args.nucs, args.make_tab_delimited, args.no_collapse = 0, True, True

    # d_dict is of the form sequence -> barcode -> record id
    d_dict = defaultdict(defaultdict)
    idx = SeqIO.index(args.infile[0], args.format) # for O(1) access to the file
    
    # Iterate through the records and put them in d_dict
    for record_id in idx:
        # get the record
        record = idx[record_id]
        
        # reset nucs (the number of nucleotides to use as the barcode)
        nucs = args.nucs

        # If we're trimming from the three-prime end, just reverse the sequence
        # (and switch back at the end)
        if args.three_prime:
            record = record[::-1]

        # get the sequence: use upper case for everything
        seq = str(record.seq).upper()
        
        # Deal with cytosines added by the RT at the 3' end of the 
        # DNA (5' end of the transcript)
        # 
        # If we require guanines, and there aren't guanines, continue
        if seq[nucs:nucs+args.require_guanines] != 'G'*args.require_guanines:
            continue
                            
        # if trimming cytosines and/or Ns as well, use them to count cDNAs
        nucs_to_remove = []
        if args.remove_guanines: nucs_to_remove.append('G')
        if args.remove_ns: nucs_to_remove.append('N')
        try: 
            while seq[nucs] in nucs_to_remove:
                nucs += 1
        except IndexError:
            if args.debug:
                sys.stderr.write('nucs:' + str(nucs) + '\n')
                sys.stderr.write('seq:' + str(seq) + '\n')
                sys.stderr.write('nucs_to_remove:' + str(nucs_to_remove) + '\n')
                raise
            else:
                continue # usually a zero-length FASTQ entry, skip it

        # split the sequence
        barcode, sequence = seq[0:nucs], seq[nucs:]

        # If --ns-are-gs, replace Ns with Gs in the barcode sequence only
        # This is useful if the sequencer is running into a ton of Gs 
        # all at the same time, and the whole flow cell lights up
        if args.ns_are_gs:
            barcode = barcode.replace('N', 'G')
        
        # ignore barcodes with Ns if no_ns is specified
        if args.no_ns and 'N' in barcode:
            continue

        # if no_collapse, we're not trimming, so each barcode gets a unique ID
        # this requires the record ids to be unique for each sequence
        if args.no_collapse:
            barcode = (barcode, record_id)
            
        # if the barcode isn't it the dictionary, or if
        # the new sequence has better quality than the old one, use it
        if barcode not in d_dict[sequence] or mean_quality(record[nucs:]) > \
            mean_quality(idx[d_dict[sequence][barcode]][nucs:]):
            d_dict[sequence][barcode] = record.id
            
    # Iterate through all the sequences and start writing them to files 
    for sequence, inner_dict in d_dict.iteritems():

        # keep track of the best record for the given barcode
        best_record, best_qual, best_barcode = None, 0, None

        # iterate through all the barcodes and records
        for barcode, record_id in inner_dict.iteritems():
            if type(barcode) is tuple:
                barcode = barcode[0]
            
            # Get the record
            record = idx[record_id]

            # If we used --three-prime, trim the three-prime end of the sequence
            if args.three_prime:
                record = record[:-len(barcode)]
            else:
                record = record[len(barcode):]
            
            # write to barcode file if we're writing all the barcodes
            # or if it matches the specific sequence
            if args.barcode_file is not None and (
                    args.barcode_sequence is None or \
                    args.barcode_sequence.upper() in sequence):
                barcode_rec = SeqRecord(Seq(barcode, IUPACAmbiguousDNA()))
                args.barcode_file.write(barcode_rec.format('fasta'))

            # If we're somehow aggregating the data, we need to find the best record to output
            if args.use_counts or args.make_tab_delimited:
                # find the best record, trim it
                record = idx[record_id][len(barcode):]
                qual = mean_quality(record)
                if qual > best_qual:
                    best_record, best_qual, best_barcode = record, qual, barcode

            # Otherwise, just write it now
            else:
                if args.append_guanines_to_header: 
                    import flulib
                    record = flulib.add_metadata_to_header(record, 
                                                           {'5ptrimmed': barcode[args.nucs:]})
                args.outfile.write(record.format(args.format))
        num_reads = len(inner_dict)
        if args.make_tab_delimited:
            args.outfile.write(str(best_record.seq) + '\t' + \
                               str(num_reads) + '\n')
        if args.use_counts:
            import flulib
            best_record = flulib.add_metadata_to_header(best_record, {'count': str(num_reads)})
            if args.append_guanines_to_header:
                best_record = flulib.add_metadata_to_header(best_record, 
                                                            {'5ptrimmed': barcode[args.nucs:]})
            best_record.name, best_record.description = record.id, record.id
            args.outfile.write(best_record.format(args.format))

def main(argv):
    'Count transcripts.'
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=\
                                     argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--no-ns', action='store_true', help="If --no-ns, " + \
                        "don't use Ns in the barcode sequence and discard " + \
                        "reads with barcodes that contain Ns.")
    parser.add_argument('-d', '--debug', action='store_true')
    parser.add_argument('infile', nargs=1, type=str,
                        help='Fastq or fasta file.')
    parser.add_argument('-o', '--outfile', nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout, help='Output file.')
    parser.add_argument('-n', '--nucs', nargs='?', type=int, default=8,
                        help='Number of nucleotides that comprise the barcode.')
    parser.add_argument('-g', '--remove-guanines', action='store_true', default=False,
                        help='If --remove-guanines, removes guanines ' + \
                        'as well and uses them for cDNA counting.')
    parser.add_argument('-u', '--append-guanines-to-header', action='store_true',
                        default=False, help='If --append-guanines-to-header, ' + \
                        'adds a count of the number of guanines to the ' + \
                        'header.')
    parser.add_argument('-x', '--remove-ns', action='store_true', default=False,
                        help='If --remove-ns, removes unknown Ns in addition ' + \
                        'to the barcode and uses them for cDNA counting.')
    parser.add_argument('-r', '--require-guanines', type=int,
                        default=0, help='If --require-guanines, will ' + \
                        'discard reads that do not contain at least n guanines' + \
                        'at the 5\' end.')
    parser.add_argument('-c', '--use-counts', action='store_true',
                        help='If --use-counts, collapse all reads and add ' + \
                        ' counts= to the read ID.')
    parser.add_argument('--barcode-file', type=argparse.FileType('w'),
                        nargs='?', default=None, help='If --barcode-file is' + \
                        ' specified, all barcodes will be written in fasta ' + \
                        'format to this file.')
    parser.add_argument('--barcode-sequence', default=None, nargs='?',
                        type=str, help='If --barcode-sequence is specified, ' + \
                        'only barcodes that match the particular downstream ' + \
                        'sequence will be written to barcode file. For ' + \
                        'example, if the barcode is ACGT and the downstream ' + \
                        'sequence is AGGTAC, if --barcode-sequence AGGTA, ' + \
                        'then write ACGT to --barcode-file.')
    parser.add_argument('--make-tab-delimited', action='store_true',
                        help='If using this option, outputs a ' + \
                        'tab-delimited file with just two columns: '+ \
                        'sequence and count number.')
    parser.add_argument('--ns-are-gs', action='store_true',
                        help='If using this option, converts Ns in the ' + \
                        'barcode to Gs (since one of the sequencing runs had ' + \
                        'lots of Gs')
    parser.add_argument('-f', '--format', nargs='?', type=str,
                        default='fastq-illumina',
                        help='''One of the file formats from Biopython SeqIO.
    Described further here: http://biopython.org/wiki/SeqIO#File_Formats''')
    parser.add_argument('--no-collapse', action='store_true',
                        help="Don't collapse reads, just trim them.")
    parser.add_argument('--make-table', action='store_true',
                        help='Equivalent: --make-tab-delimited ' + \
                        '--no-collapse -n 0. Basically make a R-compatible ' + \
                        'fastq file. ')
    parser.add_argument('-3', '--three-prime', action='store_true', default=False,
                        help='Enabling --three-prime causes count_transcripts to'
                        'collapse on the left instead of the right.')
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == "__main__":
    main(sys.argv[1:])
