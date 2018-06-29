''''Takes as input a sequence file, and then randomly shuffles all of the sequences.
Several options available to modify the extent of shuffling. 
'''

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet
import random
import itertools
import argparse
import math
from collections import defaultdict

def random_permutation(iterable, r=None):
    "Random selection from itertools.permutations(iterable, r)"
    pool = tuple(iterable)
    r = len(pool) if r is None else r
    return tuple(random.sample(pool, r))

def num_possible_permutations(iterable):
    '''If there are n items and r of them are nondistinguishable, the number
    of arrangements of these n items will be n! / (r_1! r_2! ... r_n!)
    '''
    ans = math.factorial(len(iterable)) # n!
    d = defaultdict(int)
    for elt in iterable:
        d[elt] += 1
    for r in d.itervalues():
        ans /= math.factorial(r)
    return ans

def different_num_cpgs(record1, record2):
    return str(record1.seq).count('CG') != str(record2.seq).count('CG')    

def record_from_indices(record, indices):
    '''Given a list of integers and a sequence record, will create a new
    record corresponding to the indices specified.'''
    new_seq = Seq(''.join([record[i] for i in indices]), 
                  SingleLetterAlphabet())
    new_record = SeqRecord(new_seq, record.id)
    new_record.description, new_record.name = record.id, record.id
    if record.letter_annotations:
        new_annotations = [record.letter_annotations['phred_quality'][i]
                               for i in indices]
        new_record.letter_annotations['phred_quality'] = new_annotations
    return new_record

def passes_tests(new_record, old_record, args):
    'Should we accept the shuffled sequence?'
    return not (args.preserve_cpgs and different_num_cpgs(old_record, new_record) \
        or args.not_same and str(new_record.seq) == str(old_record.seq) or \
        args.no_gs and str(new_record.seq).upper().startswith('G'))

def process_record(original_record, args):
    length = len(original_record)
    indices = range(length)
    seq = str(original_record.seq)
    # maximum possible number of permutations
    max_possible = num_possible_permutations(seq)

    # Initialize variables
    s = set()
    permutation = random_permutation(indices)
    s.add(permutation)
    record = record_from_indices(original_record, permutation)

    while not passes_tests(record, original_record, args):
        # If we've exhausted all the possibilities, return and go to next record
        if len(s) >= max_possible:
            return
        # Otherwise, create a new permutation and record, and continue
        else:
            permutation = random_permutation(indices)
            if permutation not in s:
                record = record_from_indices(original_record, permutation)
                s.add(permutation)

    # If we've gotten here, the while loop has exited, so it passed the tests
    # and we can print the record
    args.outfile.write(record.format(args.format))
    
        
def process_file(args):
    for record in SeqIO.parse(args.infile, args.format):
        process_record(record, args)

def main(argv):
    parser = argparse.ArgumentParser(argv, description=__doc__,
                                     formatter_class=\
                                     argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'), 
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'), 
                        default=sys.stdout)
    parser.add_argument('--preserve-cpgs', default=False, action='store_true')
    parser.add_argument('--not-same', default=False, action='store_true',
                        help='If --not-same, the shuffled sequence will not '
                        'be the same as the original.')
    parser.add_argument('-g', '--no-gs', default=False, action='store_true')
    parser.add_argument('-f', '--format', type=str, default='fastq-illumina')
    args = parser.parse_args(argv)
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
