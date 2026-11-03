"""SI figure (fig:SI_event_time_auc_methods): yearly per-drug probability by method, held-out (test) pool.

Per decision year, the mean over drugs of the probability that a method scores the drug's new trial disease above a
matched disease not yet trialed with the drug (ties one half; random = 0.5), with 95% CIs from 2,000 drug-bootstrap
replicates (data/derived/fig7_prospective/auc_per_drug.tsv, yearly scopes). Panel a: popularity- and
connectivity-matched alternatives (the Fig 7 comparison). Panel b: popularity-matched only. Right-hand labels give the
pooled per drug-year value. Shared-gene, GO biological-process and HPO phenotype similarity are drawn dashed with open
markers because one current annotation release is used for every decision year; DisAtlas, weighted Adamic-Adar and
popularity are rebuilt yearly. Plotted years: 2000 to 2021, matching Fig 7b.
Usage: python figures/SI/si_event_time_auc_methods.py [OUTDIR]
"""
import os, sys, json
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import sys as _s, pathlib as _p; _s.path.insert(0, str(_p.Path(__file__).resolve().parents[2] / "src")); import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE; os.makedirs(OUT, exist_ok=True)
D = os.path.join(HERE, '..', '..', 'data', 'derived', 'fig7_prospective'); START = 2000
mpl.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 6, 'axes.linewidth': .6, 'xtick.major.width': .6,
                     'ytick.major.width': .6, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
CD = '#1F4E99'; GR = '#555555'
A = pd.read_csv(f'{D}/auc_per_drug.tsv', sep='\t'); A['scope'] = A.scope.astype(str); A = A[A['pool'] == 'test']
# method, label, colour, rebuilt yearly
METH = [('disatlas_ens5', 'DisAtlas', CD, True),
        ('gene_jaccard', 'Shared genes*', '#D55E00', False),
        ('adamic_adar_w', 'Weighted Adamic–Adar', GR, True),
        ('popularity', 'Popularity alone', '#A0A0A0', True),
        ('gobp_simgic', 'GO process*', '#009E73', False),
        ('hpo_simgic', 'HPO phenotype*', '#CC79A7', False)]
COLS = [('pct_p1', 'Popularity- and connectivity-matched alternatives'), ('pct_s3', 'Popularity-matched alternatives')]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.6), sharey=True)
vals = {}
for ax, (an, title) in zip(axs, COLS):
    ends = []
    for k, (m, lab, col, yearly) in enumerate(METH):
        Z = A[(A.analysis == an) & (A.method == m) & A.scope.str.isdigit()].copy(); Z['yr'] = Z.scope.astype(int)
        Z = Z[Z.yr >= START].sort_values('yr'); assert len(Z) == 2021 - START + 1 - (START <= 1998), (an, m, len(Z))
        P = A[(A.analysis == an) & (A.method == m) & (A.scope == 'pooled_drug_year')]; assert len(P) == 1, (an, m)
        pv, plo, phi = float(P.auc.iloc[0]), float(P.lo.iloc[0]), float(P.hi.iloc[0]); vals[f'{an}|{m}'] = (pv, plo, phi, int(P.n_drugs.iloc[0]))
        z = 6 if m == 'disatlas_ens5' else 3; lw = 1.0 if m == 'disatlas_ens5' else .7
        yy = np.arange(Z.yr.min(), Z.yr.max() + 1); vv = Z.set_index('yr').auc.reindex(yy).values
        ax.vlines(Z.yr, Z.lo, Z.hi, color=col, lw=.6, alpha=.8, zorder=z)
        ax.plot(yy, vv, color=col, lw=lw, ls='-' if yearly else (0, (3, 1.5)), zorder=z)
        ax.plot(Z.yr, Z.auc, 'o', color=col, ms=2.4 if m != 'disatlas_ens5' else 2.8, mfc=col if yearly else 'white', mew=.7, zorder=z + .5)
        ends.append([pv, f'{pv:.3f}  {lab}', col])
    ax.hlines(0.5, START - .4, 2021.4, color='k', lw=.7, ls=(0, (4, 2)), zorder=1)
    ends.append([0.5, '0.500  random', 'k'])
    ends.sort(key=lambda e: e[0]); pos = []
    for e in ends:                                                     # de-overlap right-hand labels, by line count
        h = 0.019 * (e[1].count('\n') + 1)
        pos.append(e[0] if not pos else max(e[0], pos[-1][0] + pos[-1][1] / 2 + h / 2 + 0.004)); pos[-1] = (pos[-1], h)
    for (y, h), (v, t, col) in zip(pos, ends):
        ax.text(2022.9, y, t, fontsize=5.5, color=col, va='center', ha='left', linespacing=.95)
        ax.plot([2021.9, 2022.6], [v, y], color=col, lw=.5, clip_on=False)            # leader from the pooled value
    ax.text(2022.9, 0.79, 'Pooled over drug-years', fontsize=5.5, color=GR, va='top', ha='left')
    ax.set_title(title, fontsize=6.5, loc='left', pad=4)
    ax.set_xlim(START - 1, 2021.6); ax.set_xticks(range(START, 2022, 4))
    ax.set_xlabel('Decision year (evidence up to the year before)', fontsize=6)
    for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
axs[0].set_ylim(0.4, 0.8); axs[0].set_yticks([0.4, 0.5, 0.6, 0.7, 0.8])
axs[0].set_ylabel('Probability (random = 0.5)', fontsize=6)
axs[0].text(START - .6, 0.795, '\u2191 better', fontsize=5.5, color=GR, va='top', ha='left')
for ax, l in zip(axs, 'ab'): ax.text(-0.13, 1.06, l, transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom')
fig.text(0.07, 0.015, '* Current annotation release used for every decision year (not rebuilt from earlier evidence). '
         'Leaders at the right mark the pooled value.', fontsize=5.5, color=GR, ha='left', va='bottom')
fig.subplots_adjust(left=0.07, right=0.78, wspace=0.75, bottom=0.2, top=0.88)
# render-then-verify: no visible text overlaps
fig.canvas.draw(); r = fig.canvas.get_renderer()
T = [t for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
bb = [t.get_window_extent(r) for t in T]
ov = [(T[i].get_text(), T[j].get_text()) for i in range(len(T)) for j in range(i + 1, len(T)) if bb[i].overlaps(bb[j])]
assert not ov, ov
AX = [a.get_window_extent(r) for a in axs]
for t, b in zip(T, bb):
    if t.axes is axs[0]: assert not b.overlaps(AX[1]), ('crosses panel b', t.get_text())
    assert b.x1 <= fig.bbox.x1 + 1, ('beyond figure edge', t.get_text())
json.dump({'values': vals, 'start': START}, open(f'{OUT}/SI_event_time_auc_methods_values.json', 'w'), indent=1)
for ext in ('png', 'pdf'): fig.savefig(f'{OUT}/SI_event_time_auc_methods.{ext}', dpi=300, facecolor='white')
print(json.dumps({k: round(v[0], 4) for k, v in vals.items()}))
