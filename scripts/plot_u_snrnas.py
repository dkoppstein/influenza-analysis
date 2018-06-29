import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from mpltools import style
import prettyplotlib as ppl
from compare_data import SEQUENCE_COLORMAP

OUTFILE = '/lab/solexa_bartel/koppstein/influenza/new/graphs/U_snRNA_vs_literature.png'

LITERATURE_VALUES = [
    ('U1', 1000000),
    ('U2', 500000),
    ('U3', 200000),
    ('U4', 200000),
    ('U5', 200000),
    ('U6', 400000),
    ('U7', 5000),
    ('U8', 40000),
    ('U11', 10000),
    ('U12', 5000),
    ('U13', 10000)
    ]

U_SNRNA_SEQUENCES = {
    'ATCGCTTCTCG': 'U2', # U2
    'ATCGCTTCTC': 'U2', # U2
    'ATACTTACCT': 'U1', # U1
    'ATACTTACCTG': 'U1', # # U1
    'AAGACTATATTTT': 'U3', # U3
    'AAGACTATACTTTCA': 'U3', #U3
    'ATCGTCAGGT': 'U8', # U8
    'ATCGTCAGGTG': 'U8', # U8
    'ATCCTTTTGTA': 'U13', # U13
    'AGCTTTGC': 'U4', # U4 # 'AGCTTTGCGCA'
    'ATACTCTGGTTTCT': 'U5',
    'ATACTCTGGTTT': 'U5', 
    'AGTGTTACA': 'U7',
    'AAAAAGGGCTTCT': 'U11',
    'AAAAAAGGGCTTCT': 'U11',
    'ATGCCTTAAACTTA': 'U12',
    'ATGCCTTAAACTTAT': 'U12',
}
        
INFILE = '/lab/solexa_bartel/koppstein/influenza/newest/intermediate/NS1_TS_par_min3_max9_trim4_for_analysis_summary.txt'

df = pd.read_table(INFILE, index_col=None)
df.columns = [c.upper() for c in df.columns]

d = defaultdict(int)
for seq, snrna in U_SNRNA_SEQUENCES.iteritems():
    df2 = df[df['SEQUENCE'] == seq]
    d[snrna] += int(df2['COUNT'].sum())
    
counts = pd.Series(d, name='counts')
index = [x[0] for x in LITERATURE_VALUES]
literature = pd.Series(dict(LITERATURE_VALUES), name='literature')
both = pd.DataFrame(literature).reset_index().merge(pd.DataFrame(counts).reset_index(), how='outer').set_index('index')
both = both.fillna(1)
both = both.reindex(index)

style.use('ggplot')

fig = plt.figure()
ax = plt.gca()

ax.set_yscale('symlog')
ax.set_xscale('symlog')
plt.xlabel('Counts')
plt.ylabel('Literature values (copies per cell)')
plt.xlim([-1, 10000000])
plt.ylim([100, 10000000])

snrnas = [x[0] for x in LITERATURE_VALUES]

for snrna in snrnas:
    for seq, snrna_two in U_SNRNA_SEQUENCES.iteritems():
        if snrna == snrna_two:
            break
    if snrna != 'U6':
        color = SEQUENCE_COLORMAP[seq]
    else:
        color = '#b15928'
    plt.scatter(both.ix[snrna, 'counts'], both.ix[snrna, 'literature'],
                facecolor=color, edgecolor=color,
                label=snrna)

l = plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                ncol=3, borderaxespad=0., fancybox=True,
                shadow=True, scatterpoints=1, fontsize=13)

for text in l.get_texts():
    text.set_color('black')
try:
    bbox_extra_artists = (l,)
except NameError:
    bbox_extra_artists = None

plt.savefig(OUTFILE, format='png', dpi=1200,
            bbox_extra_artists=bbox_extra_artists, bbox_inches='tight')
