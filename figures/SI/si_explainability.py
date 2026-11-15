"""Supplementary Figure S26 (SI_explainability): evidence behind DisAtlas and TxGNN predictions for three rare-common pairs.

The TxGNN columns of the panel a table and the panel b sub-networks are the stored composite of the TxGNN Explorer
session (data/derived/si_explainability/SI_explainability_base.png; TxGNN published checkpoint, no model is run here).
This script redraws only the DisAtlas score column of the panel a table (header cell and three case cells, located from
the table's ruling lines) and the panel b column heading, at a font size matched to the measured ink height of the text it
replaces; every other pixel is kept.
Scores: cosine of the released gene-anchored ensemble (mean over 5 seeds), percentile among all 84,766,710 disease pairs
(upper triangle, 13,021 diseases) and rank of the common disease in the rare disease's row (1 = highest).
Input: data/derived/si_explainability/case_scores_cosine.tsv. When the released matrix
data/zenodo/model_outputs/cosine_matrix/cosine_ens5_13021.npy (with ids_13021.json) is present the three rows are
recomputed from it and must match.

    python figures/SI/si_explainability.py
Writes SI_explainability.png and SI_explainability_values.json next to this script.
"""
import os, sys, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
D = os.path.join(REPO, 'data', 'derived', 'si_explainability')
SRC = os.path.join(D, 'SI_explainability_base.png')
SC = os.path.join(REPO, 'data', 'zenodo', 'model_outputs', 'cosine_matrix')
X = pd.read_csv(os.path.join(D, 'case_scores_cosine.tsv'), sep='\t').set_index('tag')
CASES = ['Marfan-TAA', 'MODY1-T2D', 'GaucherSapC-PD24']

def recompute():
    """Recompute the three case rows from the released cosine matrix (optional; 0.7 GB)."""
    m = os.path.join(SC, 'cosine_ens5_13021.npy')
    if not os.path.exists(m): return False
    ids = json.load(open(os.path.join(SC, 'ids_13021.json'))); ix = {d: i for i, d in enumerate(ids)}; n = len(ids)
    C = np.load(m, mmap_mode='r'); allc = np.sort(np.asarray(C, np.float64)[np.triu_indices(n, 1)])
    for t in CASES:
        r = X.loc[t]; i, j = ix[r.id_a], ix[r.id_b]; c = float(C[i, j]); rc = np.asarray(C[i], np.float64).copy(); rc[i] = -np.inf
        assert round(c, 4) == r.cosine and round(100 * np.searchsorted(allc, c, 'right') / len(allc), 4) == r.cosine_pct_all_pairs
        assert int((rc > rc[j]).sum() + 1) == r.rank_b_in_row_a and len(allc) == r.n_pairs, t
    return True

im = Image.open(SRC).convert('RGB'); W, H = im.size; A = np.asarray(im.convert('L'))
dark = A < 110
# ruling lines of the panel a table (top 20% of the composite)
yt = int(0.20 * H)
rows = np.where(dark[:yt, int(0.3 * W):int(0.7 * W)].mean(1) > 0.9)[0]
hl = [g[0] for g in np.split(rows, np.where(np.diff(rows) > 2)[0] + 1)]
assert len(hl) >= 5, hl
y_top, y_bot = hl[0], hl[-1]
cols = np.where(dark[y_top:y_bot, :].mean(0) > 0.9)[0]
vl = [int(g[0]) for g in np.split(cols, np.where(np.diff(cols) > 2)[0] + 1)]
assert len(vl) >= 6 and 250 < vl[4] - vl[3] < 400, vl   # left, |case|, |rare|, |related|, |DisAtlas (then TxGNN lines)
x0, x1 = vl[3], vl[4]
row_lines = [y for y in hl if dark[y, x0 + 3:x1 - 3].mean() > 0.9]   # top, header bottom, 3 row bottoms
assert len(row_lines) == 5, (row_lines, hl)

def ink_bbox(y0, y1, xa, xb):
    s = dark[y0:y1, xa:xb]; r = np.where(s.any(1))[0]; c = np.where(s.any(0))[0]
    return (y0 + r[0], y0 + r[-1], xa + c[0], xa + c[-1]) if len(r) else None

def line_bands(y0, y1, xa, xb):
    r = np.where(dark[y0:y1, xa:xb].any(1))[0]; return [(y0 + g[0], y0 + g[-1]) for g in np.split(r, np.where(np.diff(r) > 1)[0] + 1)]

def render(text, w, h, fs, ha, va, xy):
    fig = plt.figure(figsize=(w / 300, h / 300), dpi=300); ax = fig.add_axes([0, 0, 1, 1]); ax.axis('off'); ax.set_xlim(0, w); ax.set_ylim(h, 0)
    ax.text(*xy, text, fontsize=fs, ha=ha, va=va, linespacing=1.18); fig.canvas.draw()
    a = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy(); plt.close(fig); return a

def ink_h_of(text, fs):
    a = render(text, 1200, 300, fs, 'left', 'top', (5, 5)); r = np.where((a.min(2) < 110).any(1))[0]; return r[-1] - r[0] + 1

pad = 4
b1 = line_bands(row_lines[1] + pad, row_lines[2] - pad, x0 + pad, x1 - pad)   # font size from the first line of case 1
h_digit = b1[0][1] - b1[0][0] + 1
fs = 10 * h_digit / ink_h_of('0.950', 10)
vals = [dict(tag=t, cosine=round(float(X.loc[t, 'cosine']), 3), pct=float(X.loc[t, 'cosine_pct_all_pairs']), rank=int(X.loc[t, 'rank_b_in_row_a'])) for t in CASES]
new = np.asarray(im).copy()
hy0, hy1 = row_lines[0] + pad, row_lines[1] - pad + 1   # header cell
w, h = x1 - x0 - 2 * pad, hy1 - hy0
new[hy0:hy1, x0 + pad:x1 - pad] = render('DisAtlas\ncosine\nsimilarity\n(percentile, rank)', w, h, fs, 'center', 'center', (w / 2, h / 2))
for k, v in enumerate(vals):
    cy0, cy1 = row_lines[1 + k] + pad, row_lines[2 + k] - pad + 1
    bb = ink_bbox(cy0, cy1, x0 + pad, x1 - pad); assert bb, k
    w, h = x1 - x0 - 2 * pad, cy1 - cy0
    t = f"{v['cosine']:.3f}\n({v['pct']:.2f}$^{{\\mathrm{{th}}}}$, {v['rank']})"
    a = render(t, w, h, fs, 'left', 'top', (0, 0)); r = np.where((a.min(2) < 110).any(1))[0]; c = np.where((a.min(2) < 110).any(0))[0]
    dy, dx = (bb[0] - cy0) - r[0], (bb[2] - x0 - pad) - c[0]   # anchor where the replaced text started
    new[cy0:cy1, x0 + pad:x1 - pad] = render(t, w, h, fs, 'left', 'top', (dx, dy))
# panel b column heading over the DisAtlas graphs, left quarter, first ink band below the "b" title row
yb = hl[-1] + 10; xa, xb = int(0.03 * W), int(0.26 * W)
g = [bd for bd in line_bands(yb, int(yb + 0.1 * H), xa, xb) if bd[1] - bd[0] > 8]
hb0, hb1 = g[1] if len(g) > 1 else g[0]
ob = ink_bbox(hb0, hb1 + 1, xa, xb); cx = (ob[2] + ob[3]) / 2
fsb = 10 * (hb1 - hb0 + 1) / ink_h_of('DisAtlas UPNA-DDAL', 10)   # ink height of the heading the base carries
by0, by1 = hb0 - 12, hb1 + 13; w, h = xb - xa, by1 - by0
new[by0:by1, xa:xb] = render('DisAtlas', w, h, fsb, 'center', 'center', (cx - xa, h / 2))
Image.fromarray(new).save(os.path.join(HERE, 'SI_explainability.png'), dpi=im.info.get('dpi', (300, 300)))
J = dict(recomputed_from_matrix=recompute(), n_diseases=int(X.n_diseases.iloc[0]), n_pairs_all=int(X.n_pairs.iloc[0]),
         cases={X.loc[t, 'case']: dict(tag=t, score=v['cosine'], pct=round(v['pct'], 1), pct_4dp=v['pct'], rank=v['rank']) for t, v in zip(CASES, vals)})
json.dump(J, open(os.path.join(HERE, 'SI_explainability_values.json'), 'w'), indent=1)
print(json.dumps(J))
