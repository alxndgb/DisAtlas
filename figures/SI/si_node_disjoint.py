"""Supplementary Figure SI_node_disjoint: DisAtlas with benchmark diseases seen in training vs node-disjoint at every step.

Appendix section: Node-disjoint evaluation (#sec:SI_node_disjoint).
  a  AUROC (95% disease-cluster CI) by group of IxIDN-Matched pairs for both training regimes;
  b  paired change, node-disjoint minus seen (open = CI includes 0).
Input:  figures/SI/source_data/node_disjoint_auroc.tsv
Output: figures/SI/node_disjoint_robustness.{png,pdf}
Run:    python figures/SI/si_node_disjoint.py
"""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import save, src
T = pd.read_csv(src('node_disjoint_auroc.tsv'), sep='\t')
SC = [('All 4364 pairs', 'All pairs'), ('One endpoint unseen', 'One disease unseen'), ('Both endpoints unseen', 'Both unseen'),
      ('Non-oncology pairs (label index)', 'Non-oncology'), ('Oncology pairs (label index)', 'Oncology'),
      ('Held-out test-test pairs', 'Test split'), ('Held-out test-test, both unseen', 'Test split, both unseen')]
C = {'transductive': '#4C72B0', 'node_disjoint': '#DD8452'}
NM = {'transductive': 'Benchmark diseases seen in training', 'node_disjoint': 'Node-disjoint at every step'}
plt.rcParams.update({'font.size': 8, 'axes.spines.top': False, 'axes.spines.right': False})
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw=dict(width_ratios=[1.3, 1], wspace=0.08), sharey=True)
y = np.arange(len(SC))[::-1]
for j, k in enumerate(['transductive', 'node_disjoint']):
    r = T[T.key == k].set_index('scope').loc[[s for s, _ in SC]]
    a1.errorbar(r.auroc, y + (0.15 if j == 0 else -0.15), xerr=[r.auroc - r.ci_lo, r.ci_hi - r.auroc], fmt='o', ms=4,
                color=C[k], ecolor=C[k], elinewidth=1, capsize=0, label=NM[k])
a1.set_yticks(y); a1.set_yticklabels([n for _, n in SC]); a1.set_xlabel('AUROC (95% disease-cluster CI)')
a1.axvline(0.5, color='0.7', lw=0.6, ls=':'); a1.set_xlim(0.58, 0.92)
fig.legend(*a1.get_legend_handles_labels(), frameon=False, fontsize=7, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1.0)); a1.set_title('a  AUROC by group of pairs', loc='left', fontsize=8, fontweight='bold')
r = T[T.key == 'node_disjoint'].set_index('scope').loc[[s for s, _ in SC]]
sig = (r.d_hi < 0) | (r.d_lo > 0)
a2.hlines(y, r.d_lo, r.d_hi, color=C['node_disjoint'], lw=1)
a2.scatter(r.delta_vs_transductive, y, s=22, color=np.where(sig, C['node_disjoint'], 'white'), edgecolor=C['node_disjoint'], zorder=3)
a2.axvline(0, color='0.3', lw=0.7); a2.set_xlabel('Node-disjoint minus seen (paired 95% CI)')
a2.set_title('b  Paired change', loc='left', fontsize=8, fontweight='bold'); a2.margins(x=0.08)
print(save(fig, 'node_disjoint_robustness', bbox_inches='tight'))

# ---- value check: appendix #sec:SI_node_disjoint (3 d.p.) ----
_v = T.set_index(['key', 'scope'])
_r = lambda k, s, c='auroc': round(float(_v.loc[(k, s), c]) + 1e-12, 3)
assert (_r('node_disjoint', 'All 4364 pairs'), _r('transductive', 'All 4364 pairs'), _r('node_disjoint', 'All 4364 pairs', 'delta_vs_transductive')) == (0.749, 0.765, -0.016)
assert (_r('node_disjoint', 'Both endpoints unseen'), _r('transductive', 'Both endpoints unseen'), _r('node_disjoint', 'Both endpoints unseen', 'delta_vs_transductive')) == (0.733, 0.766, -0.033)
assert (_r('node_disjoint', 'Held-out test-test pairs'), _r('transductive', 'Held-out test-test pairs')) == (0.758, 0.783)
print('appendix values: OK')
