'''
This works like fastx_collapser, except that it additionally requires metadata 
present in the header to be the same for sequences that it collapses. 

Outputs a fasta file. 

For example, 

>arbitrarystringANNOTATION:5ptrimmed=G|3ptrimmed=GCA
ACGTACGTACGT

>foobarANNOTATION:5ptrimmed=G|3ptrimmed=GCA
ACGTACGTACGT

Will get collapsed to 

>foobarANNOTATION:5ptrimmed=G|3ptrimmed=GCA|count=2
ACGTACGTACGT

but if another sequence is present, but has different metadata, 
say, 

>foobarANNOTATION:5ptrimmed=G|3ptrimmed=GCAGGG
ACGTACGTACGT

it will *not* be collapsed. 
'''

import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import flulib
import argparse
from collections import defaultdict

def process_file(args):
    # (sequence, header_set) -> count
    seq_dict = defaultdict(int)
    for record in SeqIO.parse(args.infile, args.format):
        # Create a hashable, frozen set of tuples of the header dict
        # 
        # We require a unique sequence and header set for each record 
        # to be counted
        header_set = frozenset(tuple((k.upper(),v) \
        for k,v in flulib.headerdict_from_string(record.id).iteritems()))
        seq = str(record.seq)
        seq_dict[seq, header_set] += 1
    num_seq = 0
    for ((seq, header_set), count) in seq_dict.iteritems():
        # Make a new dictionary from the header and add the count attribute
        new_dict = dict(header_set)
        new_dict.update({'COUNT': str(count)})
        record = flulib.add_metadata_to_header(
            SeqRecord(Seq(seq), id=str(num_seq) + '-'), new_dict)
        # We can't keep qualities, so make a fasta file
        args.outfile.write(record.format('fasta'))
        num_seq += 1

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-f', '--format', type=str, nargs='?',
                        default='fastq-illumina')
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
