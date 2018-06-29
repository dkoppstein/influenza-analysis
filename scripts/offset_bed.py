import pybedtools
import argparse
import sys

'''Offsets the intervals of BED files in a strand-specific manner. 
Accepts positive negative numbers on the command line.''' 

def process_file(args):
    for iv in pybedtools.BedTool(args.infile):
        if args.offset_three_prime_wrt_five_prime == 0:
            if iv.strand == '+':
                new_start = int(iv.start) + args.five_prime_offset
                new_stop = int(iv.stop) + args.three_prime_offset
            else:
                new_stop = int(iv.stop) - args.five_prime_offset
                new_start = int(iv.start) - args.three_prime_offset
            # If this creates entries below zero, throw them out
        else:
            new_start = int(iv.start)
            if iv.strand == '+':
                new_stop = int(iv.start) + \
                  args.offset_three_prime_wrt_five_prime + 1
            else:
                new_stop = int(iv.stop)
                new_start = int(iv.stop) - \
                    args.offset_three_prime_wrt_five_prime - 1
        if new_start < 1: continue
        if new_stop < 1: continue
        iv.start, iv.stop = new_start, new_stop
        args.outfile.write(str(iv))

def main(argv):
    parser = argparse.ArgumentParser(argv)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'), 
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', default=sys.stdout, 
                        type=argparse.FileType('w'))
    parser.add_argument('--five-prime-offset', type=int, nargs='?', default=0)
    parser.add_argument('--three-prime-offset', type=int, nargs='?', default=0)
    parser.add_argument('--offset-three-prime-wrt-five-prime', 
                        type=int, nargs='?', default=0)
    args = parser.parse_args()
    process_file(args)

if __name__ == '__main__':
    main(sys.argv[1:])
