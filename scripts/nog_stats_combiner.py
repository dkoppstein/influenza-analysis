
import pandas as pd
import sys
import argparse

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infiles', type=str, nargs='+')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    args = parser.parse_args(argv)
    gs, nogs = 0, 0
    for f in args.infiles:
        df = pd.read_table(f, index_col=None)
        gs += df['Gs'][0]
        nogs += df['No_Gs'][0]
    args.outfile.write('Gs\tNo_Gs\n')
    args.outfile.write('%s\t%s\n' % (gs, nogs))

if __name__ == '__main__':
    main(sys.argv[1:])
