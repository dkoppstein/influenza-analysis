import argparse
import sys
import pandas as pd
from cStringIO import StringIO
from subprocess import Popen, PIPE

def intersect_tables(table1, table2):
    'Returns a tuple of tables containing only the reads represented by both.'
    table1_set, table2_set = (set(table['READ_NAME']) for table in
                              (table1, table2))
    intersect = table1_set.intersection(table2_set)
    return tuple(restrict_table_to_set(table, intersect) for table in 
                 (table1, table2))

def restrict_table_to_set(table, s):
    'Only returns a row from table if in set s'
    return pd.DataFrame((row for __, row in table.iterrows() if
                        row['READ_NAME'] in s), columns=table.columns)

def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-a', '--file-a', type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-b', '--file-b', type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('-oa', '--outfile-a', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-ob', '--outfile-b', type=argparse.FileType('w'),
                        required=True)
    parser.add_argument('--gzipped', action='store_true')
    args = parser.parse_args(argv)
    if args.gzipped:
        args.file_a, args.file_b = (StringIO(Popen(['zcat', f.name],
                                          stdout=PIPE).communicate()[0]) \
                                          for f in [args.file_a, args.file_b])
    table1, table2 = tuple((pd.read_table(f, index_col=None) \
                           for f in (args.file_a, args.file_b)))
    table1, table2 = intersect_tables(table1, table2)
    for table, outfile in [(table1, args.outfile_a),
                            (table2, args.outfile_b)]:
        table.to_csv(outfile, sep='\t', index=False, na_rep='NA')
    
if __name__ == '__main__':
    main(sys.argv[1:])
