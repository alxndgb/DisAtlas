"""Supplementary Figure SI_ablation_aggregation: aggregation-rule ablation.

Appendix section: Aggregation-rule ablation (#sec:SI_ablations_aggregation).
  a  share of DisAtlas network edges kept by each alternative rule vs its AUROC on all 4,364 IxIDN-Matched pairs (95% CI);
  b  paired change in AUROC vs the mean-rank rule on all, test-split and non-oncology pairs (open = CI includes 0).
Inputs (figures/SI/source_data/): ablation_retrained_arms_auroc.tsv, aggregation_edges_kept.json
Output: figures/SI/ablation_aggregation_rule.{png,pdf}
Run:    python figures/SI/si_ablation_aggregation.py
"""
import os, json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import save, src
T = pd.read_csv(src('ablation_retrained_arms_auroc.tsv'), sep='\t'); T = T[T.metric == 'auroc']
O = json.load(open(src('aggregation_edges_kept.json')))['arms']
SC = [('IxIDN-Matched (all 4,364 pairs)', 'All pairs', 'o'), ('Test-split pairs', 'Test-split', 'D'), ('Non-oncology pairs', 'Non-oncology', 's')]
A = [('control', 'Mean rank\n(DisAtlas)'), ('crank', 'CRank'), ('mean', 'Mean similarity'), ('max', 'Max similarity'), ('rra', 'Kolde RRA')]
COL = {'control': '#222222', 'crank': '#d17a22', 'mean': '#1b6ca8', 'max': '#2a9d8f', 'rra': '#8d6cab'}
g = lambda k, s: T[(T.key == k) & (T.scope == s)].iloc[0]
fig, (a, b) = plt.subplots(1, 2, figsize=(8.2, 3.4), gridspec_kw=dict(width_ratios=[1, 1.35], wspace=0.35))
c = g('control', SC[0][0]); a.axhspan(c.lo, c.hi, color='0.9', zorder=0)
for k, nm in A:
    r = g(k, SC[0][0]); x = 100.0 if k == 'control' else O[k]['pct_control_edges_kept']
    a.errorbar(x, r.value, yerr=[[r.value - r.lo], [r.hi - r.value]], fmt='o', color=COL[k], ms=6, capsize=2, lw=1)
    dx, dy, ha = {'control': (-7, 0, 'right'), 'mean': (-7, 0, 'right'), 'max': (-7, 9, 'right'), 'crank': (7, -9, 'left'), 'rra': (7, 0, 'left')}[k]
    a.annotate(nm, (x, r.value), xytext=(dx, dy), textcoords='offset points', fontsize=7, color=COL[k], ha=ha, va='center', bbox=dict(fc='white', ec='none', pad=0.6, alpha=0.85), zorder=5)
a.set_xlim(40, 106); a.set_xlabel('Edges kept vs. DisAtlas (%)', fontsize=8); a.set_ylabel('AUROC, all 4,364 pairs (95% CI)', fontsize=8)
a.tick_params(labelsize=7.5)
ys = np.arange(4)[::-1]
for j, (s, lab, mk) in enumerate(SC):
    off = (1 - j) * 0.22
    for i, (k, nm) in enumerate(A[1:]):
        r = g(k, s); sig = r.d_lo > 0 or r.d_hi < 0
        b.plot([r.d_lo, r.d_hi], [ys[i] + off] * 2, color=COL[k], lw=1.1, alpha=1 - 0.25 * j)
        b.plot(r.d_control, ys[i] + off, mk, ms=5, color=COL[k] if sig else 'white', mec=COL[k], mew=1.1)
b.axvline(0, color='0.3', lw=0.8); b.set_yticks(ys); b.set_yticklabels([nm for _, nm in A[1:]], fontsize=7.5)
b.set_xlabel('Change in AUROC vs mean rank (95% CI)', fontsize=8); b.tick_params(axis='x', labelsize=7.5)
from matplotlib.lines import Line2D
b.legend(handles=[Line2D([], [], marker=mk, ls='', color='0.3', label=lab) for _, lab, mk in SC] +
         [Line2D([], [], marker='o', ls='', mfc='white', mec='0.3', label='CI includes 0')], fontsize=6.5, loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=4, frameon=False)
for ax in (a, b):
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
a.set_title('a', loc='left', fontweight='bold', fontsize=10); b.set_title('b', loc='left', fontweight='bold', fontsize=10)
print(save(fig, 'ablation_aggregation_rule', bbox_inches='tight'))

# ---- value check: appendix #sec:SI_ablations_aggregation (3 d.p.) ----
_r = lambda x: round(float(x) + 1e-12, 3)
assert _r(g('control', SC[0][0]).value) == 0.765
r = g('crank', SC[0][0]); assert (_r(r.value), _r(r.d_control), _r(r.d_lo), _r(r.d_hi)) == (0.747, -0.017, -0.033, -0.003)
assert _r(g('rra', SC[0][0]).d_control) == -0.013
assert sorted(round(O[k]['pct_control_edges_kept']) for k in O)[::3] == [61, 79]
for k in ('crank', 'mean', 'max', 'rra'):
    r = g(k, SC[1][0]); assert abs(_r(r.d_control)) <= 0.018 and r.d_lo < 0 < r.d_hi, k
print('appendix values: OK')
