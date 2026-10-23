"""Supplementary Figure SI_ablation_multimetric: multi-metric component ablation on the test split.

Appendix section: Component ablations (#fig:SI_ablation_multimetric). Rows: the DisAtlas model, the retrained
combination-rule and annotation-source variants and three conventional similarity baselines; columns: AUROC, AUPRC,
precision in the top 5% / 10%, per-disease Hits@1, F1 at the top-half cut, then benchmark and network coverage.
Cell colour = change against the retrained rank-mean control; bold = paired 95% disease-cluster bootstrap CI
(B = 1,000) excludes zero. Values are the stored results of the retrained variants; the DisAtlas row is recomputed
here from data/08S_eval/ixidn4364_scores.tsv.
Inputs: source_data/ablation_multimetric_test_split.tsv, source_data/ablation_multimetric_coverage.tsv,
  data/08S_eval/ixidn4364_scores.tsv
Output: figures/SI/ablation_multimetric_test.{png,pdf}
Run:    python figures/SI/si_ablation_multimetric.py
"""
import os, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import save, src, REPO_ROOT

T = pd.read_csv(src('ablation_multimetric_test_split.tsv'), sep='\t').set_index(['key', 'metric'])
COV = pd.read_csv(src('ablation_multimetric_coverage.tsv'), sep='\t').set_index('key')
for a, b, p in [('pairs_scored', 'n_pairs', 'pct_pairs'), ('bench_diseases_scored', 'n_bench_diseases', 'pct_bench_diseases'),
                ('net_diseases', 'n_net', 'pct_net')]:
    COV[p] = (100 * COV[a] / COV[b]).round(2)
MET = ['auroc', 'auprc', 'p_top5', 'p_top10', 'hits1', 'f1bal']
XB, XD, BW = len(MET) - 0.35, len(MET) + 1.15, 1.35   # left edge of the two coverage columns, column width
MLAB = ['AUROC', 'AUPRC', 'Precision\ntop 5%', 'Precision\ntop 10%', 'Hits@1\nper disease', 'F1\n(top-half cut)']
GROUPS = [('DisAtlas', ['frozen']),
          ('Combination rule (retrained)', ['control', 'crank', 'mean', 'max', 'rra']),
          ('Annotation sources (retrained)', ['go_hpo', 'go_only', 'hpo_only']),
          ('Conventional similarity, no network', ['hpo_go_resnik_fusion', 'go_resnik_genebma', 'mondo_simgic'])]
LABEL = {'frozen': 'DisAtlas', 'control': 'Rank-mean, 4 sources (control)', 'crank': 'CRank', 'mean': 'Mean of similarities',
         'max': 'Max of similarities', 'rra': 'Kolde RRA', 'go_hpo': 'GO + HPO only', 'go_only': 'GO only', 'hpo_only': 'HPO only',
         'hpo_go_resnik_fusion': 'HPO+GO Resnik fusion', 'go_resnik_genebma': 'GO Resnik gene-BMA', 'mondo_simgic': 'MONDO SimGIC'}


def prec_top(s, y, f):
    o = np.argsort(-s, kind='stable'); k = np.arange(1, len(s) + 1) <= f * len(s) + 1e-9
    return y[o][k].sum() / k.sum()


def disatlas_test_metrics():
    """Point estimates of the DisAtlas row and the caption counts, from the shipped pair scores."""
    P = pd.read_csv(os.path.join(REPO_ROOT, 'data', '08S_eval', 'ixidn4364_scores.tsv'), sep='\t')
    D = P[P.pool == 'test'].reset_index(drop=True); y = D.y.values.astype(float); s = D.p.values
    q = pd.concat([pd.DataFrame({'d': D.a, 'i': D.index}), pd.DataFrame({'d': D.b, 'i': D.index})], ignore_index=True)
    q['y'] = y[q.i]; q['s'] = s[q.i]
    g = q.groupby('d').y.agg(['sum', 'count']); elig = g.index[(g['sum'] > 0) & (g['sum'] < g['count'])]
    e = q[q.d.isin(elig)]; e = e[e.s == e.groupby('d').s.transform('max')]
    o = np.argsort(-s, kind='stable'); pred = np.zeros(len(s)); pred[o[np.arange(1, len(s) + 1) <= 0.5 * len(s) + 1e-9]] = 1
    m = dict(auroc=roc_auc_score(y, s), auprc=average_precision_score(y, s), p_top5=prec_top(s, y, .05),
             p_top10=prec_top(s, y, .10), hits1=e.groupby('d').y.mean().mean(), f1bal=f1_score(y, pred))
    return len(D), len(set(D.a) | set(D.b)), len(elig), m


n_pairs, n_dis, n_hits1, m = disatlas_test_metrics()
assert (n_pairs, n_dis, n_hits1) == (1094, 447, 228), (n_pairs, n_dis, n_hits1)
for k, v in m.items():
    assert abs(v - T.loc[('frozen', k), 'value']) < 5e-6, (k, v, T.loc[('frozen', k), 'value'])
assert (COV.n_pairs == 1094).all() and (COV.n_bench_diseases == 447).all() and (COV.n_net == 13021).all()
assert set(T.index.get_level_values(0)) == {k for _, g in GROUPS for k in g}

norm = TwoSlopeNorm(vmin=-0.25, vcenter=0.0, vmax=0.05); cmap = plt.get_cmap('RdBu')
nrow = sum(len(g) + 1 for _, g in GROUPS); H = 0.26 * (nrow + 3) + 0.6
fig, ax = plt.subplots(1, 1, figsize=(11.0, H), gridspec_kw=dict(hspace=0.11))
fig.subplots_adjust(bottom=0.13)
yy = 0
for gname, keys in GROUPS:
    ax.text(-0.5, yy, gname, fontsize=8.5, fontweight='bold', va='center', ha='left'); yy += 1
    for k in keys:
        ax.text(-0.6, yy, LABEL[k], ha='right', va='center', fontsize=8, fontweight='bold' if k == 'frozen' else 'normal')
        for j, mt in enumerate(MET):
            r = T.loc[(k, mt)]; sig = (r.d_lo > 0 or r.d_hi < 0) and k not in ('control', 'frozen')
            ax.add_patch(plt.Rectangle((j - 0.47, yy - 0.42), 0.94, 0.84, color=cmap(norm(np.clip(r.d_control, -0.25, 0.05))), lw=0))
            ax.text(j, yy, f'{r.value:.3f}', ha='center', va='center', fontsize=7.5, color='white' if r.d_control < -0.15 else 'black',
                    fontweight='bold' if sig else 'normal')
        cv = COV.loc[k]; bc = '#1b4f8a' if k == 'frozen' else '#4a7fb5'
        for x0, pct, txt in ((XB, cv.pct_pairs, f'{cv.pct_pairs:.1f}% ({cv.pct_bench_diseases:.1f}%)'),
                             (XD, cv.pct_net, f'{cv.pct_net:.2f}%' if 99.9 < cv.pct_net < 100 else f'{cv.pct_net:.1f}%')):
            ax.add_patch(plt.Rectangle((x0, yy - 0.36), BW, 0.72, color=(0.94, 0.94, 0.94), lw=0))
            ax.add_patch(plt.Rectangle((x0, yy - 0.36), BW * pct / 100, 0.72, color=bc, lw=0))
            ax.text(x0 + 0.05 if pct > 55 else x0 + BW * pct / 100 + 0.04, yy, txt, ha='left', va='center', fontsize=6.6,
                    color='white' if pct > 55 else 'black')
        if k == 'frozen':
            ax.add_patch(plt.Rectangle((-0.5, yy - 0.5), len(MET), 1.0, fill=False, ec='black', lw=1.0))
        yy += 1
    yy += 0.4
ax.set_xlim(-0.5, XD + BW + 0.05); ax.set_ylim(yy - 0.3, -1.6); ax.set_yticks([])
ax.set_xticks(list(range(len(MET))) + [XB + BW / 2, XD + BW / 2])
ax.set_xticklabels(MLAB + ['Benchmark coverage\n% pairs scored\n(% diseases)', 'DisAtlas coverage\n% of 13,021\nnetwork diseases'], fontsize=8)
ax.xaxis.tick_top(); ax.tick_params(length=0)
for s_ in ax.spines.values(): s_.set_visible(False)
ax.set_title('Test-split pairs (1,094)', loc='left', fontsize=10, fontweight='bold', x=-0.21, y=1.07)
cax = fig.add_axes([0.32, 0.06, 0.35, 0.02])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation='horizontal', ticks=[-0.2, -0.1, 0, 0.05])
cb.ax.tick_params(labelsize=7); cb.set_label('Change vs control (bold: paired 95% disease-cluster CI excludes 0)', fontsize=7.5)
save(fig, 'ablation_multimetric_test', bbox_inches='tight'); plt.close(fig)
