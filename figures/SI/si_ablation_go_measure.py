"""Supplementary Figure SI_ablation_gomeasure: gene-pair GO distance rule.

Appendix section: GO gene-pair distance rule (#sec:SI_ablations_gomeasure).
  a  closest-decile gene pairs under the minimum rule that fall above the median distance under the mean over term pairs
     or the best-match average, per GO namespace;
  b  minimum-rule network links kept per alternative rule, at the matched cutoff and at cutoff 0.2
     (in-bar labels: links kept | IxIDN-Matched pairs scored | network diseases placed);
  c  change in AUROC vs the minimum rule (95% CI): matched cutoff, cutoff 0.2, and cutoff 0.2 on common-scorable pairs.
Inputs (figures/SI/source_data/): go_measure_auroc.tsv, go_measure_common_scorable_auroc.tsv, go_measure_network_topology.tsv,
  go_measure_rule_disagreement.json, go_measure_diseases_placed.tsv
Output: figures/SI/ablation_go_measure.{png,pdf}
Run:    python figures/SI/si_ablation_go_measure.py
"""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import save, src
rd = src
T = pd.read_csv(rd('go_measure_auroc.tsv'), sep='\t'); T = T[T.scope == 'All 4364 pairs'].set_index('key')
C = pd.read_csv(rd('go_measure_common_scorable_auroc.tsv'), sep='\t'); C = C[C.scope.str.startswith('All')].set_index('key')
P = pd.read_csv(rd('go_measure_network_topology.tsv'), sep='\t').set_index('key')
D = json.load(open(rd('go_measure_rule_disagreement.json')))['namespaces']
M5 = ['lca_bma', 'lca_mean', 'resnik_max', 'resnik_bma', 'lin_bma']
ML = {'lca_bma': 'LCA BMA', 'lca_mean': 'LCA mean', 'resnik_max': 'Resnik max', 'resnik_bma': 'Resnik BMA', 'lin_bma': 'Lin BMA'}
BL, RE, LR = '#1b6ca8', '#b2182b', '#e8a09a'
fig, ax = plt.subplots(1, 3, figsize=(11.5, 3.4), gridspec_kw=dict(width_ratios=[0.8, 1.3, 1.25], wspace=0.5))
a = ax[0]; ns = [('biological_process', 'BP'), ('molecular_function', 'MF'), ('cellular_component', 'CC')]; x = np.arange(3)
a.bar(x - 0.18, [100 * D[n]['frac_close_by_min_but_mean_above_median'] for n, _ in ns], 0.36, color='0.6', label='Mean over term pairs')
a.bar(x + 0.18, [100 * D[n]['frac_close_by_min_but_bma_above_median'] for n, _ in ns], 0.36, color='#2a9d8f', label='Best-match average')
a.set_xticks(x); a.set_xticklabels([l for _, l in ns], fontsize=8); a.set_ylabel('Closest-decile gene pairs (min)\nabove median distance (%)', fontsize=8)
a.legend(fontsize=6.5, frameon=False, loc='upper right')
b = ax[1]; y = np.arange(len(M5))[::-1]
b.barh(y + 0.19, [P.loc['qmatch_' + m, 'pct_kept'] for m in M5], 0.38, color=BL, label='Matched cutoff')
b.barh(y - 0.19, [P.loc['literal_' + m, 'pct_kept'] for m in M5], 0.38, color=RE, label='Cutoff 0.2')
# coverage per arm (percent): IxIDN-Matched pairs scored (4,364 minus unembedded) and network diseases placed
# (network diseases the retrained model embeds, i.e. with >=1 training edge: go_measure_diseases_placed.tsv)
EMB = pd.read_csv(rd('go_measure_diseases_placed.tsv'), sep='\t').set_index('arm')
fmtp = lambda v: f'{v:.2f}%' if 99.9 < v < 100 else f'{v:.1f}%'
for off, pre in ((0.19, 'qmatch_'), (-0.19, 'literal_')):
    for i, m in enumerate(M5):
        k = pre + m; pp = 100 * (4364 - T.loc[k, 'n_unembedded']) / 4364; pd_ = 100 * EMB.loc[k, 'embedded_in_13021'] / 13021
        b.text(1.5, y[i] + off, f"{P.loc[k, 'pct_kept']:.1f}%  |  {pp:.1f}%  |  {fmtp(pd_)}", fontsize=6.3, va='center', ha='left', color='white')
b.set_yticks(y); b.set_yticklabels([ML[m] for m in M5], fontsize=8); b.set_xlim(0, 101); b.axvline(100, color='0.3', lw=0.7, ls=':')
b.set_xticks([0, 25, 50, 75, 100])
b.set_xlabel('Minimum-rule links kept (%)\nin bars: links kept | IxIDN-Matched pairs scored | network diseases placed', fontsize=7.5); b.legend(fontsize=6.5, frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
c = ax[2]
for off, src, pre, col, lab in [(0.25, T, 'qmatch_', BL, 'Matched cutoff'), (0, T, 'literal_', RE, 'Cutoff 0.2'), (-0.25, C, 'literal_', LR, 'Cutoff 0.2, common-scorable pairs')]:
    for i, m in enumerate(M5):
        r = src.loc[pre + m]; c.plot([r.d_lo, r.d_hi], [y[i] + off] * 2, color=col, lw=1.2); c.plot(r.delta_vs_control, y[i] + off, 'o', color=col, ms=4, label=lab if i == 0 else None)
c.axvline(0, color='0.3', lw=0.8); c.set_yticks(y); c.set_yticklabels([ML[m] for m in M5], fontsize=8)
c.set_xlabel('Change in AUROC vs minimum rule (95% CI)', fontsize=8); c.legend(fontsize=6.5, frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.2), ncol=2)
for i, A in enumerate(ax):
    A.tick_params(labelsize=7.5); A.set_title('abc'[i], loc='left', fontweight='bold', fontsize=10)
    for sp in ('top', 'right'): A.spines[sp].set_visible(False)
print(save(fig, 'ablation_go_measure', bbox_inches='tight'))

# ---- value check: appendix #sec:SI_ablations_gomeasure ----
_r = lambda x: round(float(x) + 1e-12, 3)
assert (round(100 * D['biological_process']['frac_close_by_min_but_mean_above_median'], 1),
        round(100 * D['biological_process']['frac_close_by_min_but_bma_above_median'], 1)) == (22.7, 1.5)
q = ['qmatch_' + m for m in M5]
assert (P.loc[q, 'pct_kept'].round(1).min(), P.loc[q, 'pct_kept'].round(1).max()) == (99.7, 99.8)
assert (T.loc[q, 'delta_vs_control'].round(3).min(), T.loc[q, 'delta_vs_control'].round(3).max()) == (-0.009, -0.003)
assert _r(T.loc['lca_min', 'auroc']) == 0.764
# fixed cutoff 0.2: the appendix ranges (65.7-70.8% links kept, AUROC -0.209 to -0.174, 37.8-62.8% pairs scored)
# cover the four measures other than LCA BMA (94.0% kept, -0.017, 99.6% scored; see panels b-c)
l4 = ['literal_' + m for m in M5 if m != 'lca_bma']
pp = 100 * (4364 - T.loc[l4, 'n_unembedded']) / 4364
assert (P.loc[l4, 'pct_kept'].round(1).min(), P.loc[l4, 'pct_kept'].round(1).max()) == (65.7, 70.8)
assert (T.loc[l4, 'delta_vs_control'].round(3).min(), T.loc[l4, 'delta_vs_control'].round(3).max()) == (-0.209, -0.174)
assert (pp.round(1).min(), pp.round(1).max()) == (37.8, 62.8)
print('appendix values: OK')
