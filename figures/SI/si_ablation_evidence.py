"""Supplementary Figure SI_ablation_evidence: one-component-at-a-time ablation of DisAtlas on IxIDN-Matched.

Appendix sections: Component Ablations (#sec:SI_ablations), Evidence-source ablation (#sec:SI_ablations_evidence),
Graph encoder (#sec:SI_ablations_encoder).
  a  paired change in AUROC vs DisAtlas (95% CI) for evidence-source, representation, decoder and aggregation arms on all
     4,364 pairs (unscored pairs ranked last), on the pairs each arm can score (open circles) and on the 1,094 test-split pairs.
  b  benchmark coverage (% of pairs scored, % of 938 diseases) and DisAtlas coverage (% of 13,021 network diseases).
Inputs (figures/SI/source_data/): ablation_retrained_arms_auroc.tsv, ablation_direct_distance_baselines.csv,
  ablation_decoder_vs_cosine.csv, ablation_arm_graph_coverage.csv, ablation_coverage.tsv
Output: figures/SI/ablation_evidence_sources.{png,pdf}
Run:    python figures/SI/si_ablation_evidence.py
"""
import os, pandas as pd, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import save, src
T = pd.read_csv(src('ablation_retrained_arms_auroc.tsv'), sep='\t'); T = T[T.metric == 'auroc']
DB = pd.read_csv(src('ablation_direct_distance_baselines.csv')); DV = pd.read_csv(src('ablation_decoder_vs_cosine.csv'))
CV = pd.read_csv(src('ablation_arm_graph_coverage.csv')).set_index('arm')
# disease coverage of the direct distances: same annotation requirement and scorable pairs as the matching
# conventional measure in ablation_coverage.tsv (asserted below on pair counts)
COV = pd.read_csv(src('ablation_coverage.tsv'), sep='\t').set_index(['scope', 'key']).loc['IxIDN-Matched (all 4,364 pairs)']
DMATCH = {'mechdist_go': 'go_resnik_genebma', 'mechdist_hpo': 'hpo_resnik_bma'}
SC = {'all': 'IxIDN-Matched (all 4,364 pairs)', 'test': 'Test-split pairs'}
def arm(k, s):
    r = T[(T.key == k) & (T.scope == SC[s])].iloc[0]
    sp = (r.auroc_scorable - r.control_auroc_same_pairs) if r.pairs_scored < 4364 else np.nan
    return r.d_control, r.d_lo, r.d_hi, sp
def direct(col, s):
    r = DB[(DB.column == col) & (DB.scope == {'all': 'all', 'test': 'test-test'}[s])].iloc[0]
    return -r.delta_all, -r.dall_ci_hi, -r.dall_ci_lo, -r.delta_model_minus
def cosine(s):
    r = DV[(DV.iloc[:, 0] == {'all': 'all', 'test': 'test-test'}[s]) & (DV.iloc[:, 1].astype(str) == 'ens') & (DV.iloc[:, 2] == 'cosine')].iloc[0]
    d, lo, hi = r.iloc[-3], r.iloc[-2], r.iloc[-1]   # decoder - cosine
    return -d, -hi, -lo, np.nan
ROWS = [('Evidence sources', [('GO + HPO (no HPA, Pfam)', 'go_hpo'), ('GO + HPA + Pfam (no HPO)', 'gene_only'), ('HPA + Pfam', 'px_only'),
                              ('GO only', 'go_only'), ('HPA only', 'hpa_only'), ('Pfam only', 'dom_only'), ('HPO only', 'hpo_only')]),
        ('Representation', [('No network: direct GO distance†', ('direct', 'mechdist_go')), ('No network: direct HPO distance†', ('direct', 'mechdist_hpo')),
                            ('GCN autoencoder', 'gcn'), ('GraphSAGE autoencoder', 'sage'), ('Laplacian eigenmap', 'lap')]),
        ('Decoder', [('Embedding cosine, no decoder†', ('cos',))]),
        ('Aggregation', [('CRank (submitted rule)', 'crank'), ('Mean of similarities', 'mean'), ('Max of similarities', 'max'), ('Kolde RRA', 'rra')])]
lab, grp, V = [], [], {'all': [], 'test': []}
for g, items in ROWS:
    for nm, k in items:
        lab.append(nm); grp.append(g)
        for s in V:
            V[s].append(direct(k[1], s) if isinstance(k, tuple) and k[0] == 'direct' else cosine(s) if isinstance(k, tuple) else arm(k, s))
COL = {'Evidence sources': '#1b6ca8', 'Representation': '#2a9d8f', 'Decoder': '#8d6cab', 'Aggregation': '#d17a22'}
n = len(lab); _y, cur, prev = [], 0.0, None
for g_ in grp:
    if prev is not None and g_ != prev: cur += 0.8
    _y.append(cur); cur += 1; prev = g_
y = -np.array(_y)
fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 6.4), gridspec_kw=dict(width_ratios=[2.3, 1.15], wspace=0.04), sharey=True)
for s, off, mk in (('all', 0.14, 'o'), ('test', -0.14, 'D')):
    for i, (d, lo, hi, sp) in enumerate(V[s]):
        c = COL[grp[i]]
        a.plot([lo, hi], [y[i] + off] * 2, color=c, lw=1.2, alpha=0.9 if s == 'all' else 0.55)
        a.plot(d, y[i] + off, mk, color=c, ms=5 if s == 'all' else 4, alpha=1 if s == 'all' else 0.6, mec='white', mew=0.4)
        if s == 'all' and not np.isnan(sp): a.plot(sp, y[i] + off, 'o', mfc='white', mec=c, ms=5, mew=1.1)
a.axvline(0, color='0.3', lw=0.8); a.set_yticks(y); a.set_yticklabels(lab, fontsize=7.5)
a.set_xlabel('Change in AUROC vs DisAtlas (95% CI)', fontsize=8); a.tick_params(axis='x', labelsize=7.5)
a.grid(axis='x', color='0.92', lw=0.6); a.set_axisbelow(True)
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([], [], marker='o', ls='', color='0.3', label='All 4,364 pairs (unscored pairs last)'),
                  Line2D([], [], marker='o', ls='', mfc='white', mec='0.3', label='Pairs the variant can score'),
                  Line2D([], [], marker='D', ls='', color='0.3', alpha=0.6, ms=4, label='Test-split pairs (1,094)')],
         fontsize=6.5, loc='lower left', frameon=False)
yy = 0
for g, items in ROWS:
    top = y[yy]; bot = y[yy + len(items) - 1]
    a.text(1.0, (top + bot) / 2, '', transform=a.get_yaxis_transform())
    if yy: a.axhline(top + 0.9, color='0.8', lw=0.6)
    a.text(-0.335, top + 0.55, g, fontsize=7.5, fontweight='bold', color=COL[g])
    yy += len(items)
a.set_xlim(-0.34, 0.07)
# panel b: benchmark and DisAtlas coverage per arm, read from ablation_coverage.tsv
XB, XD, BW = 0.0, 1.5, 1.35
keys = [k for g, it in ROWS for _, k in it]
for i, k in enumerate(keys):
    if isinstance(k, tuple) and k[0] == 'direct':
        cv = COV.loc[DMATCH[k[1]]]; r = DB[(DB.column == k[1]) & (DB.scope == 'all')].iloc[0]
        assert int(cv.pairs_scored) == int(r.n_scorable), (k, cv.pairs_scored, r.n_scorable)
    elif isinstance(k, tuple):   # embedding cosine on the DisAtlas model
        cv = COV.loc['frozen']
    else:
        cv = COV.loc[k]
        if k in CV.index:   # trained-model coverage can only be lower than graph-presence coverage
            assert 100 - cv.pct_pairs >= CV.loc[k, 'pct_pairs_lost'] - 0.05, k
    for x0, pct, txt in ((XB, cv.pct_pairs, f'{cv.pct_pairs:.1f}% ({cv.pct_bench_diseases:.1f}%)'),
                         (XD, cv.pct_net, f'{cv.pct_net:.2f}%' if 99.9 < cv.pct_net < 100 else f'{cv.pct_net:.1f}%')):
        b.add_patch(plt.Rectangle((x0, y[i] - 0.34), BW, 0.68, color='0.94', lw=0))
        b.add_patch(plt.Rectangle((x0, y[i] - 0.34), BW * pct / 100, 0.68, color=COL[grp[i]], lw=0))
        b.text(x0 + 0.05 if pct > 55 else x0 + BW * pct / 100 + 0.04, y[i], txt, ha='left', va='center', fontsize=6.6,
               color='white' if pct > 55 else 'black')
b.set_xlim(-0.05, XD + BW + 0.05); b.set_xticks([]); b.spines['left'].set_visible(False); b.spines['bottom'].set_visible(False)
b.tick_params(axis='y', length=0)
for x0, t in ((XB, 'Benchmark coverage\n% of 4,364 pairs scored\n(% of 938 diseases)'), (XD, 'DisAtlas coverage\n% of 13,021\nnetwork diseases')):
    b.text(x0 + BW / 2, y[0] + 0.75, t, ha='center', va='bottom', fontsize=7)
for ax in (a, b):
    for sp_ in ('top', 'right'): ax.spines[sp_].set_visible(False)
a.set_title('a', loc='left', fontsize=10, fontweight='bold', pad=30); b.set_title('b', loc='left', fontsize=10, fontweight='bold', pad=30)
print(save(fig, 'ablation_evidence_sources', bbox_inches='tight'))

# ---- value check: appendix #sec:SI_ablations_evidence (3 d.p.) ----
_r = lambda x: round(float(x) + 1e-12, 3)
_a = lambda k, s='IxIDN-Matched (all 4,364 pairs)': T[(T.key == k) & (T.scope == s)].iloc[0]
_c = COV
for k, pp, pn, au, te, sc, ct in (('go_only', 81.7, 69.0, 0.685, 0.704, 0.739, 0.768), ('hpo_only', 34.8, 78.8, 0.539, 0.569, 0.661, 0.746)):
    r = _a(k)
    assert (round(_c.loc[k, 'pct_pairs'], 1), round(_c.loc[k, 'pct_net'], 1)) == (pp, pn), k
    assert (_r(r.value), _r(_a(k, 'Test-split pairs').value), _r(r.auroc_scorable), _r(r.control_auroc_same_pairs)) == (au, te, sc, ct), k
# signal loss on own scorable pairs: appendix 0.029 (GO; difference of the rounded AUROCs, unrounded 0.028) and 0.085 (HPO)
assert round(_r(_a('go_only').control_auroc_same_pairs) - _r(_a('go_only').auroc_scorable), 3) == 0.029
assert _r(_a('hpo_only').control_auroc_same_pairs - _a('hpo_only').auroc_scorable) == 0.085
r = _a('go_hpo'); assert (_r(r.value), _r(r.d_control), _r(r.d_lo), _r(r.d_hi)) == (0.75, -0.015, -0.028, -0.002)
r = DB[(DB.column == 'mechdist_go') & (DB.scope == 'all')].iloc[0]   # direct GO distance, on the pairs it can score
assert (_r(r.delta_model_minus), _r(r.d_ci_lo), _r(r.d_ci_hi)) == (0.157, 0.108, 0.201)
r = DV[DV.scope == 'all'].iloc[0]   # decoder minus cosine, five-seed ensemble
assert (_r(-r.delta_decoder_minus_distance), _r(-r.d_ci_hi), _r(-r.d_ci_lo)) == (0.002, -0.003, 0.007)
print('appendix values: OK')
