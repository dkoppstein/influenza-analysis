#!/usr/bin/env python

'''
Use this script instead of zcat for opening high-throughput sequencing data. 
Briefly, trims the .tar header and the binary crap at the end of the file. 
'''


import sys
import argparse
from subprocess import Popen, PIPE

def process_file(inhandle, outhandle):
    # Open a zcat subprocess and read from the inhandle
    zcat = Popen(['zcat'], stdin=inhandle, stdout=PIPE)
    # This will eventually trim off the last line
    head = Popen(['head', '-n', '-1'], stdin=PIPE, stdout=outhandle)
    # Read in first line, but only output the first part of the FASTQ file and 
    # everything after it
    first_line = zcat.stdout.readline()
    head.stdin.write(first_line[first_line.find('@'):])
    # Feed all the other lines into the head subprocess
    for line in zcat.stdout:
        head.stdin.write(line)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', nargs='?', type=argparse.FileType('r'),
                        default=sys.stdin, help='A Solexa 1.3-1.7 fastq file')
    parser.add_argument('-o', '--outfile', nargs='?', type=argparse.FileType('w'),
                        default=sys.stdout, help='Output file.')
    args = parser.parse_args()
    if args.infile.isatty():
        parser.error('No input detected.')
    process_file(args.infile, args.outfile)

if __name__ == '__main__':
    main()
