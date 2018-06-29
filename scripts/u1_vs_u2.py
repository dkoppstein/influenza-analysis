import pandas as pd
import argparse
import sys
from plot_u_snrnas import U_SNRNA_SEQUENCES
from compare_data import SEQUENCE_COLORMAP
import numpy as np
from collections import defaultdict
from cStringIO import StringIO
from subprocess import Popen, PIPE
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

def process_files(args):
    gen_one_lst, gen_two_lst = [], []
    colormap = {}
    for lst, f_lst in [(gen_one_lst, args.gen_one_ts_infiles),
                       (gen_two_lst, args.gen_two_ts_infiles)]:
        for f in f_lst:
            if args.gzipped:
                f = StringIO(Popen(['zcat', f.name], stdout=PIPE
                                   ).communicate()[0])
            d = defaultdict(int)
            df = pd.read_table(f, index_col=None)
            df.columns = [s.upper() for s in df.columns]
            for seq, snrna in U_SNRNA_SEQUENCES.iteritems():
                if snrna == 'U1' or snrna == 'U2':
                    df2 = df[df['SEQUENCE'] == seq]
                    d[snrna] += float(df2['COUNT'].sum()) * 1000000. / float(df['COUNT'].sum())
                    colormap[snrna] = SEQUENCE_COLORMAP[seq]
            counts = pd.Series(d, name=f)
            lst.append(counts)
    gen_one_df, gen_two_df = (pd.DataFrame(l).T for l in (gen_one_lst, gen_two_lst))
    gen_one_means, gen_two_means = (df.apply(np.mean, axis=1) for df in (gen_one_df, gen_two_df))
    gen_one_stds, gen_two_stds = (df.apply(np.std, axis=1) for df in (gen_one_df, gen_two_df))
    (((gen_one_u1_mean, gen_one_u2_mean), (gen_two_u1_mean, gen_two_u2_mean)),
     ((gen_one_u1_std, gen_one_u2_std), (gen_two_u1_std, gen_two_u2_std))) = \
      (((s['U1'], s['U2']) for s in series) \
      for series in [(gen_one_means, gen_two_means), (gen_one_stds, gen_two_stds)])
    z_stat, p_val = ttest_ind(gen_one_df.ix['U1', :], gen_two_df.ix['U2', :])
    sys.stderr.write('%s%s%s%s' % (gen_one_u1_mean, gen_one_u2_mean,
                                   gen_two_u1_mean, gen_two_u2_mean))
    fig, ax = plt.subplots()

    index = np.arange(1)
    bar_width = 0.1

    opacity = 0.4

    # U1
    rects1 = ax.bar(index+bar_width, [gen_one_u1_mean], bar_width, color=colormap['U1'], yerr=[gen_one_u1_std], ecolor='black')
    # U2
    rects2 = ax.bar(index+2*bar_width, [gen_one_u2_mean], bar_width, color=colormap['U2'], yerr=[gen_one_u2_std], ecolor='black')

    ax.legend((rects1[0], rects2[0]), ('U1', 'U2'))
    plt.annotate('p = %f (t-test)' % p_val,
                 xy=(0.05, 0.95), xycoords='axes fraction')
    ax.get_xaxis().set_ticks([])
    plt.xlim(0, 0.4)
    plt.ylabel('Counts per million')
    plt.savefig(args.outfile, format=args.format)
    
def main(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gen-one-ts-infiles', nargs='+',
                        type=argparse.FileType('r'))
    parser.add_argument('--gen-two-ts-infiles', nargs='+',
                        type=argparse.FileType('r'))
    parser.add_argument('--gzipped', action='store_true')
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                        default=sys.stdout)
    parser.add_argument('-f', '--format', default='png')
    args = parser.parse_args(argv)
    process_files(args)

if __name__ == '__main__':
    main(sys.argv[1:])
