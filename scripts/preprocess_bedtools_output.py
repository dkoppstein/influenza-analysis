import argparse
import sys
import flulib
from deparse_bedtools_output import (
     BEDTOOLS_INTERSECT_COLUMNS, COLUMNS_TO_ADD, GTF_KEYS
     )

'''
This file takes the output of "bedtools intersect" and turns it into a
pandas-friendly TSV file, with uppercase column names. 
'''                                     

def deparse_gtf_metadata(metadata):
    'Returns a dict of key, value pairs from the GTF metadata string.'
    metadata = metadata.replace('"', '')
    if 'zero_length_insertion' in metadata:
        metadata = metadata[:metadata.find('zero_length_insertion')].strip()
    metadata = (tuple(i.split()) for i in metadata.split('; '))
    d = dict(metadata)
    d = {k.upper(): v for k, v in d.iteritems()}
    return d

def process_file(args):
    total_list = BEDTOOLS_INTERSECT_COLUMNS + COLUMNS_TO_ADD + GTF_KEYS
    args.outfile.write('\t'.join(total_list) + '\n')
    for line in args.infile:
        new_dict = dict(zip(total_list, ['NaN']*len(total_list)))
        line = line.strip().split('\t')
        line_dict = dict(zip(BEDTOOLS_INTERSECT_COLUMNS, line))
        new_dict.update(line_dict)
        new_dict.update(flulib.headerdict_from_string(line_dict['READ_NAME']))
        new_dict.update(deparse_gtf_metadata(line_dict['TSS_META']))
        new_dict = {k.upper() : v for k, v in new_dict.iteritems()} # uppercase
        final_line = '\t'.join([new_dict[v] for v in total_list])
        args.outfile.write(final_line + '\n')

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--infile', type=argparse.FileType('r'),
                        default=sys.stdin)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args(argv)
    process_file(args)
    
if __name__ == '__main__':
    main(sys.argv[1:])
