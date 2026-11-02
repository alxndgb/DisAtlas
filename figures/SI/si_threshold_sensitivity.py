#!/usr/bin/env python
"""SI figure: threshold sensitivity of the GO/HPO distance cutoff.

The network was rebuilt at each cutoff and retrained with three seeds; scoring IxIDN-Matched needs the
retrained embeddings, so the per-cutoff results ship as precomputed tables:
  data/derived/threshold_sensitivity/threshold_sweep_values.tsv    per-cutoff graph size and AUROC (disease-cluster
                                                                   bootstrap, B = 2000, 95% CI)
  data/derived/threshold_sensitivity/threshold_sweep_per_seed.tsv  per-seed AUROC
"""
import sys as _s, pathlib as _p; _s.path.insert(0, str(_p.Path(__file__).resolve().parents[2] / "src")); import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
import sys as _sys; _sys.path.insert(0, str(REPO / "src"))
DV = REPO / "data" / "derived" / "threshold_sensitivity"
V = pd.read_csv(DV / "threshold_sweep_values.tsv", sep="\t").sort_values("cutoff").reset_index(drop=True)
N_PAIRS = int(V.n_all.iloc[0]); assert N_PAIRS == 4364
CHOSEN = '0.20'

# ---------------------------------------------------------------- figure (layout of the earlier SI figure)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.spines.top': False, 'axes.spines.right': False,
                     'axes.linewidth': 0.8, 'xtick.major.width': 0.8, 'ytick.major.width': 0.8})
DARK, LIGHT, SHADE, CHC, GREY = '#1f4e9c', '#6a9fd8', '#f9efe3', '#d9822b', '#555555'
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.1), gridspec_kw=dict(width_ratios=[1.55, 1], wspace=0.28))
x = V.cutoff.values; i0 = int(np.argmin(np.abs(x - 0.10)))
for ax in (a1, a2):
    ax.axvspan(0.175, 0.225, color=SHADE, zorder=0, lw=0); ax.set_xticks(x); ax.set_xticklabels([f'{v:.2f}' for v in x])
    ax.set_xlabel('GO/HPO distance cutoff', fontsize=8)
a1.plot(x, V.auroc_test, '--s', color=LIGHT, lw=1.4, ms=5, zorder=3)
a1.plot(x, V.auroc_all, '-o', color=DARK, lw=1.8, ms=5, zorder=4)
if V.n_unscorable[i0] > 0:
    a1.plot(x[i0], V.auroc_all[i0], 'o', mfc='white', mec=DARK, mew=1.4, ms=5, zorder=5)
    a1.annotate(f'0.10: {V.n_unscorable[i0]:,} pairs\nunscorable', xy=(x[i0], V.auroc_all[i0]), xytext=(x[i0] + 0.006, V.auroc_all[i0] - 0.03),
                fontsize=6, color=GREY, arrowprops=dict(arrowstyle='-', color=GREY, lw=0.6))
lo = min(V.auroc_all.min(), V.auroc_test.min()); hi = max(V.auroc_all.max(), V.auroc_test.max())
a1.set_ylim(np.floor((lo - 0.035) * 100) / 100, np.ceil((hi + 0.025) * 100) / 100); a1.set_xlim(0.085, 0.355)
a1.text(0.2, a1.get_ylim()[1], f'chosen ({CHOSEN})', ha='center', va='top', color=CHC, fontsize=6.5)
a1.text(0.307, V.auroc_test.values[-1], 'test-split\npairs', color=LIGHT, va='center', fontsize=7)
a1.text(0.307, V.auroc_all.values[-1], f'all {N_PAIRS:,}\npairs', color=DARK, va='center', fontsize=7)
off = 0.006 if V.auroc_test.values[-1] > V.auroc_all.values[-1] else -0.006
a1.set_ylabel('IxIDN-Matched AUROC', fontsize=8); a1.set_title('Retrained at each cutoff (3 seeds)', loc='left', fontsize=8)
a2.plot(x, V.n_nodes / 1000, '-o', color=GREY, lw=1.6, ms=5, zorder=3)
dn = int(V.n_nodes[V.cutoff == 0.20].iloc[0] - V.n_nodes[i0])
a2.annotate(f'\u2212{dn:,} diseases', xy=(x[i0], V.n_nodes[i0] / 1000), xytext=(x[i0] + 0.035, V.n_nodes[i0] / 1000 + 0.1),
            fontsize=6, color=GREY, arrowprops=dict(arrowstyle='-', color=GREY, lw=0.6))
a2.set_xlim(0.09, 0.31); a2.set_ylim(np.floor(V.n_nodes.min() / 500) / 2 - 0.5, np.ceil(V.n_nodes.max() / 500) / 2 + 0.3)
a2.set_ylabel('Diseases in network (thousands)', fontsize=8); a2.set_title('Only 0.10 drops diseases', loc='left', fontsize=8)
for ax, l in ((a1, 'a'), (a2, 'b')):
    ax.text(-0.13 if ax is a1 else -0.16, 1.06, l, transform=ax.transAxes, fontsize=9, fontweight='bold')
for ext in ('png', 'pdf'): fig.savefig(HERE / f'si_threshold_sensitivity.{ext}', dpi=300, bbox_inches='tight')
print('saved', HERE)
