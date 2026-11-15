"""Figure 4 / evidence-trace prose numbers (main text 'Tracing a disease-disease similarity to its evidence',
'Rheumatoid arthritis and two rare inflammatory myopathies', Methods 'Evidence trace and attribution').

Reads only data/derived/fig4_evidence_trace and writes figures/main/fig4_numbers.json: one key per manuscript value,
holding the raw value and the string as printed in the manuscript.
Percentiles are of the 5-seed ensemble cosine among all 84,766,710 disease pairs; if the released ensemble cosine
matrix is present under data/zenodo it is recomputed from it, otherwise the shipped value is used.
"""
import json, os, sys
import numpy as np, pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(ROOT, 'data', 'derived', 'fig4_evidence_trace')
OUT = os.path.join(ROOT, 'figures', 'main', 'fig4_numbers.json')
MAT = os.path.join(ROOT, 'data', 'zenodo', 'model_outputs', 'cosine_matrix', 'cosine_ens5_13021.npy')
IDS = os.path.join(ROOT, 'data', 'zenodo', 'model_outputs', 'cosine_matrix', 'ids_13021.json')

V = json.load(open(os.path.join(D, 'fig4_values.json')))
I = pd.read_csv(os.path.join(D, 'trace_displayed_items.tsv'), sep='\t')
G = pd.read_csv(os.path.join(D, 'trace_evidence_type_shares.tsv'), sep='\t')
S = pd.read_csv(os.path.join(D, 'case_cosines.tsv'), sep='\t', dtype={'seed': str})
T10 = pd.read_csv(os.path.join(D, 'trace_top10_modules.tsv'), sep='\t')
B = pd.read_csv(os.path.join(D, 'm2_vs_cos_bins.tsv'), sep='\t')

recomputed = False
if os.path.exists(MAT) and os.path.exists(IDS):
    ids = json.load(open(IDS)); ix = {d: i for i, d in enumerate(ids)}; n = len(ids)
    C = np.load(MAT, mmap_mode='r'); allc = np.sort(np.asarray(C, np.float64)[np.triu_indices(n, 1)])
    for k in S.index[S.seed == 'ensemble']:
        c = float(C[ix[S.at[k, 'id_a']], ix[S.at[k, 'id_b']]]); assert abs(c - S.at[k, 'cosine']) < 1e-5
        p = 100 * np.searchsorted(allc, c, 'right') / len(allc); assert abs(p - S.at[k, 'pct_all_pairs']) < 1e-3
    recomputed = True

out = {}
def put(key, value, fmt):
    out[key] = dict(value=value, printed=fmt.format(value))

for pair in S.pair.unique():
    s = S[S.pair == pair]; seeds = s[s.seed != 'ensemble'].cosine; e = s[s.seed == 'ensemble'].iloc[0]
    assert abs(seeds.mean() - e.cosine) < 1e-5
    put(f'{pair}.cosine_mean5', float(e.cosine), '{:.3f}')
    put(f'{pair}.cosine_seed_min', float(seeds.min()), '{:.3f}')
    put(f'{pair}.cosine_seed_max', float(seeds.max()), '{:.3f}')
    put(f'{pair}.pct_all_pairs', float(e.pct_all_pairs), '{:.1f}')

put('n_ixidn_pairs', int(V['n_pairs']), '{:,}')
put('n_pairs_with_shared_neighbour', int(V['n_m2_pos']), '{:,}')
put('spearman_logM2_cosine', float(V['spearman_logM2_cos']), '{:.2f}')
for p, k in V['n_shared_neighbours'].items(): put(f'{p}.n_shared_neighbours', int(k), '{}')
put('n_diseases_idf', int(V['idf_n_diseases']), '{:,}')
put('caption.n_equal_count_bins', int((B.bin >= 0).sum()), '{}')

def item(pair, lst, it):
    r = I[(I.pair == pair) & (I.list == lst) & (I.item == it)]; assert len(r) == 1, (pair, lst, it, len(r))
    return float(r.frac_of_M2.iloc[0])
for pair, lst, it, key in [
        ('aHUS-PNH', 'module', 'REACTOME:R-HSA-977606', 'reactome_regulation_of_complement'),
        ('aHUS-PNH', 'module', 'NAMED:complement_all', 'named_complement_all'),
        ('aHUS-PNH', 'module', 'MECHNODE:M1', 'mechanism_complement_regulator_loss'),
        ('aHUS-PNH', 'gene', 'ENSEMBL:ENSG00000125730', 'gene_C3'), ('aHUS-PNH', 'gene', 'ENSEMBL:ENSG00000106804', 'gene_C5'),
        ('aHUS-PNH', 'gene', 'ENSEMBL:ENSG00000000971', 'gene_CFH'),
        ('NMOSD-PNH', 'module', 'REACTOME:R-HSA-977606', 'reactome_regulation_of_complement'),
        ('NMOSD-PNH', 'module', 'LITMECH:MECH:M5|v1', 'litmech_terminal_pathway'),
        ('NMOSD-PNH', 'module', 'LITMECH:complement_C5_MAC|C5b-9/MAC', 'litmech_C5b9_MAC'),
        ('NMOSD-PNH', 'module', 'REACTOME:R-HSA-6798695', 'reactome_neutrophil_degranulation'),
        ('NMOSD-PNH', 'gene', 'ENSEMBL:ENSG00000125730', 'gene_C3'),
        ('Crohn-UC', 'module', 'GENEANCHOR:ENSEMBL:ENSG00000113302', 'shared_rare_gene_IL12B'),
        ('Crohn-UC', 'module', 'GENEANCHOR:ENSEMBL:ENSG00000167207', 'shared_rare_gene_NOD2')]:
    put(f'{pair}.{key}', item(pair, lst, it), '{:.3f}')

g = G[G.pair == 'aHUS-PNH'].set_index('group').share
put('aHUS-PNH.share_broad', float(g['Broad similarity']), '{:.3f}')
put('aHUS-PNH.share_shared_rare_gene', float(g['Shared rare gene']), '{:.3f}')
put('aHUS-PNH.share_modules_rest', float(1 - g['Broad similarity'] - g['Shared rare gene']), '{:.3f}')
put('RA-polymyositis.top10_items_share', float(T10[T10.pair == 'RA-polymyositis'].frac_of_M2.sum()), '{:.2f}')
put('Crohn-UC.top10_modules_share', float(T10[T10.pair == 'Crohn-UC'].frac_of_M2.sum()), '{:.2f}')
ifn = I[(I.pair == 'RA-polymyositis') & I.label.str.contains('Interferon gamma|Interleukin-10|Interleukin-4', regex=True)].frac_of_M2
put('RA-polymyositis.max_cytokine_pathway_share', float(ifn.max()), '{:.3f}')
assert ifn.max() < 0.01

json.dump(dict(recomputed_percentiles_from_matrix=recomputed, values=out), open(OUT, 'w'), indent=1)
open(OUT, 'a').write('\n')
print(len(out), 'values; percentiles recomputed:', recomputed)
