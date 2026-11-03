"""Main Figure 5: DisAtlas disease map (K = 55 clusters) and rarity-stratified benchmark performance.

Panels:
  a) Disease map (Isomap layout of DisAtlas similarities), coloured by the K = 55 clusters; six control pairs
     (four positive, two negative) joined by coloured lines; an inset shows the four eculizumab-approved diseases
     (aHUS, PNH and generalised MG in one cluster, NMOSD in the neighbouring cluster, dotted line).
  b) Common / rare disease map with five rare-common examples from different therapeutic areas.
  c) Benchmark AUROC by pair-rarity stratum: DisAtlas and TxGNN (published checkpoint; scores and protocol of the
     cross-system figure, Figure 6), one set of disease-cluster bootstrap draws.
  d) DisAtlas score distributions (true pair vs non-pair) by stratum; white dots = medians, coloured bars = control pairs.

Inputs (data/derived/):
  fig5_clustering/disease_map_k55.tsv.gz       13,021 diseases: Isomap (cosine, 15 neighbours) layout, K = 55 cluster
                                               label, disease class / rarity / name. map_x/map_y are the layout after a
                                               similarity transform (rotation/reflection, uniform scale, shift) used for
                                               display orientation only.
  fig5_clustering/ixidn_matched_pair_scores.tsv  4,364 IxIDN-Matched pairs with the DisAtlas score p.
  fig5_clustering/fig5c_txgnn_pair_scores.tsv  TxGNN score of each pair (fig5c_export_txgnn_scores.py).
  fig6_comparison/per_run_scores.tsv           disease order of the bootstrap clusters (shared with the SI).
  control_pairs/fig5_panel_a_pairs_p.tsv       DisAtlas scores of the six panel-a control pairs.
The cluster labels, layout and pair scores are regenerated from the released model by
pipeline/clustering (Zenodo pair-score matrix); this script only draws.

    python figures/main/fig5_disease_clustering_map.py [OUTDIR]     (B bootstrap draws: FIG5_B, default 2000)
Writes disease_clustering_map.png/.pdf, fig5_panel_c_auroc.tsv, fig5_panel_c_txgnn_auroc.tsv,
fig5_panel_c_paired_differences.tsv, disease_clustering_map_values.tsv.
"""
import os, sys, json, textwrap
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
DV = os.path.join(ROOT, "data", "derived")
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE

_m = pd.read_csv(os.path.join(DV, "fig5_clustering", "disease_map_k55.tsv.gz"), sep="\t")
ids = _m.mondo_id.tolist(); n = len(ids); ix = {d: i for i, d in enumerate(ids)}
assert n == 13021, n
cl = pd.DataFrame({"mds1": _m.map_x.to_numpy(), "mds2": _m.map_y.to_numpy(), "cluster": _m.cluster.to_numpy(),
                   "disease_class": _m.disease_class.to_numpy(), "is_rare": _m.is_rare.to_numpy(), "label": _m.label.to_numpy()},
                  index=pd.Index(ids, name="mondo_id"))
Y = cl[["mds1", "mds2"]].to_numpy(); clu_of = cl["cluster"]
print("n", len(cl), "K", cl.cluster.nunique())
assert cl.cluster.nunique() == 55
_c = {k: int(cl.cluster[k]) for k in ("MONDO:0009688", "MONDO:0019100", "MONDO:0100244", "MONDO:0016244")}; print("eculizumab clusters", _c)
assert _c["MONDO:0009688"] == _c["MONDO:0100244"] == _c["MONDO:0016244"] != _c["MONDO:0019100"], _c   # MG, PNH, aHUS together; NMOSD apart
_PA = pd.read_csv(os.path.join(DV, "control_pairs", "fig5_panel_a_pairs_p.tsv"), sep="\t")
PAIR_P = {(a, b): p for a, b, p in zip(_PA.a, _PA.b, _PA.p_upna_ddal)}
# ---- Panel A: six control pairs on the K = 55 clusters ----
# Eculizumab-approved diseases (PNH, aHUS, generalised MG, NMOSD) are all labelled; PNH/aHUS also form a control pair.
GOOD, BAD = "#00b140", "#b42318"
PANEL_A_PAIRS = [
    ("TSC", "MONDO:0001734", "LAM", "MONDO:0011705", "pos"),
    ("psoriasis", "MONDO:0005083", "psoriatic arthritis", "MONDO:0011849", "pos"),
    ("aHUS", "MONDO:0016244", "PNH", "MONDO:0100244", "pos"),
    ("Crohn disease", "MONDO:0005011", "ulcerative colitis", "MONDO:0005101", "pos"),
    ("diabetes insipidus", "MONDO:0004782", "T2D", "MONDO:0005148", "neg"),
    ("NF1 (mut)", "MONDO:0018208", "NF2", "MONDO:0007039", "neg"),
]
# Eculizumab inset: aHUS-PNH and PNH-gMG are joined by solid green lines (same cluster); PNH-NMOSD by a
# dotted green line (NMOSD sits in the neighbouring cluster); NMOSD is the only yellow point.
ECU_PAIRS = [("PNH", "MONDO:0100244", "gMG", "MONDO:0009688", "solid"),
             ("PNH", "MONDO:0100244", "NMOSD", "MONDO:0019100", "dotted")]
PANEL_A_SINGLES = []

# Five rare-common examples from five therapeutic areas; each rare monogenic disease sits in the cluster of a related
# common disease (asserted below).
PANEL_B_PAIRS = [  # (rare label, rare id, common label, common id)
    ("CESD (LAL-D)", "MONDO:0019149", "familial hypercholesterolemia", "MONDO:0005439"),
    ("Crigler\u2013Najjar type 2", "MONDO:0011725", "Gilbert syndrome", "MONDO:0007745"),
    ("CPVT1", "MONDO:0010975", "sudden arrhythmic death syndrome", "MONDO:0054866"),
    ("Huntington disease", "MONDO:0007739", "Parkinson disease", "MONDO:0005180"),
    ("congenital fibrosis of the\nextraocular muscles", "MONDO:0007614", "ptosis", "MONDO:0000728"),
]
for _r, _a, _c, _b in PANEL_B_PAIRS:
    assert _a in ix and _b in ix and int(cl.cluster[_a]) == int(cl.cluster[_b]), (_r, _c)
    assert cl.disease_class[_a] == "rare" and cl.disease_class[_b] == "common", (_r, _c)

import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import patheffects as pe
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import ListedColormap

from disatlas.figstyle import apply_style, panel_label   # house style
apply_style()
COMMON_COL, RARE_COL = "#2a6ebb", "#d1495b"

n_clusters = int(cl["cluster"].nunique())
cluster_ids_sorted = np.sort(cl["cluster"].unique())
palette = (list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("tab20b").colors) +
           list(plt.get_cmap("tab20c").colors))[:n_clusters]
cluster_cmap = ListedColormap(palette)
cluster_to_colidx = {cc: i for i, cc in enumerate(cluster_ids_sorted)}
point_colidx = cl["cluster"].map(cluster_to_colidx).to_numpy()


def wrap_label_capped(name, start_width=16, max_lines=2):
    w = start_width
    lines = textwrap.wrap(name, width=w, break_long_words=False)
    while len(lines) > max_lines and w < len(name):
        w += 3
        lines = textwrap.wrap(name, width=w, break_long_words=False)
    return "\n".join(lines), len(lines)


def place_leader_labels_sided(ax, items, x0, x1, y0, y1, fontsize=8.5, region_gap_frac=0.10,
                               local_pad_frac=0.09, line_height_frac=0.052, inter_item_frac=0.03,
                               push_to_edge_frac=0.0, single_region=False):
    y_mid = (y0 + y1) / 2.0 if not single_region else -np.inf
    span = y1 - y0
    for side in ("left", "right"):
        for region in ("top", "bottom"):
            sel = [it for it in items if it[3] == side and ((it[0][1] >= y_mid) == (region == "top"))]
            if not sel:
                continue
            sel = sorted(sel, key=lambda it: -it[0][1])
            n = len(sel)
            ys_actual = np.array([it[0][1] for it in sel])
            n_lines = [it[2] for it in sel]

            heights = np.array([nl * line_height_frac * span for nl in n_lines])
            gaps = inter_item_frac * span
            needed_span = heights.sum() + gaps * (n - 1)

            centre = (ys_actual.max() + ys_actual.min()) / 2.0
            edge_gap = region_gap_frac * span * 0.2
            if region == "top":
                edge_centre = (y1 - edge_gap) - needed_span / 2.0
                centre = centre + push_to_edge_frac * max(0.0, edge_centre - centre)
            else:
                edge_centre = (y0 + edge_gap) + needed_span / 2.0
                centre = centre + push_to_edge_frac * min(0.0, edge_centre - centre)

            hi_local = centre + needed_span / 2.0
            lo_local = centre - needed_span / 2.0
            if region == "top" and hi_local > y1 - edge_gap:
                shift = hi_local - (y1 - edge_gap)
                hi_local -= shift; lo_local -= shift
            if (region == "bottom" or single_region) and lo_local < y0 + edge_gap:
                shift = (y0 + edge_gap) - lo_local
                hi_local += shift; lo_local += shift

            ys = np.empty(n)
            cursor = hi_local
            for i in range(n):
                ys[i] = cursor - heights[i] / 2.0
                cursor -= heights[i] + gaps

            for (a, t, _, _), ty in zip(sel, ys):
                ax.annotate(t, a, xytext=(x0 - 0.05 if side == "left" else x1 + 0.05, ty),
                            ha="right" if side == "left" else "left", va="center", fontsize=fontsize,
                            linespacing=1.25,
                            arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3", shrinkB=4), zorder=6)


def place_labels_local(ax, items, anchors_all, fontsize=10, legend=None, avoid=()):
    """Put each label next to its own point, on the side of the map nearest that point, with a short leader.
    Greedy search over nearby offsets; a candidate is rejected if its text box overlaps a placed label, any labelled
    marker, the legend or the axes edge. items: ((x, y), text, n_lines, side)."""
    fig = ax.figure; fig.canvas.draw(); rr = fig.canvas.get_renderer(); tr = ax.transData
    abox = ax.get_window_extent(rr); taken = []
    if legend is not None: taken.append(legend.get_window_extent(rr).expanded(1.04, 1.08))
    for a_ in avoid: taken.append(a_.get_tightbbox(rr).expanded(1.04, 1.06))   # inset frame + its title
    for (x, y) in anchors_all:
        px, py = tr.transform((x, y)); taken.append(mpl.transforms.Bbox([[px - 6, py - 6], [px + 6, py + 6]]))
    for ln in ax.lines:   # connector lines: labels must not sit on them
        xd, yd = ln.get_xdata(), ln.get_ydata()
        if len(xd) != 2: continue
        for f_ in np.linspace(0, 1, 25):
            px, py = tr.transform((xd[0] + f_ * (xd[1] - xd[0]), yd[0] + f_ * (yd[1] - yd[0])))
            taken.append(mpl.transforms.Bbox([[px - 3, py - 3], [px + 3, py + 3]]))
    pxu = abs(tr.transform((1, 0))[0] - tr.transform((0, 0))[0])   # pixels per data unit
    for k in [k for k in range(len(items)) if isinstance(items[k][3], tuple)]:   # fixed-slot labels first
        (x, y), text, _, spec = items[k]; fx, fy = spec[1]
        t = ax.text(fx, fy, text, ha="left", va="center", fontsize=fontsize, linespacing=1.15, zorder=7,
                    path_effects=[pe.withStroke(linewidth=3, foreground="white")])
        bb = t.get_window_extent(rr).expanded(1.03, 1.10)
        assert abox.x0 < bb.x0 and bb.x1 < abox.x1 and abox.y0 < bb.y0 and not any(bb.overlaps(o) for o in taken), text
        if spec[0] == "slide":
            # slide the box right, 1 px at a time, until it would touch a map dot, a placed box or a connector line,
            # or come within 2.5% of the axis width of the right edge; keep the last clear position
            dots = np.vstack([c.get_offsets() for c in ax.collections if hasattr(c, "get_offsets") and len(c.get_offsets())])
            dpx = tr.transform(dots); rdot = 2.0   # dot radius in px (s=1.4 background points, plus margin)
            blk = list(taken)
            for a_ in ax.artists:   # inset connector lines (ConnectionPatch), sampled in display coords
                if hasattr(a_, "_get_path_in_displaycoord"):
                    v = a_._get_path_in_displaycoord().vertices
                    for i_ in range(len(v) - 1):
                        for f_ in np.linspace(0, 1, 40):
                            q = v[i_] + f_ * (v[i_ + 1] - v[i_]); blk.append(mpl.transforms.Bbox([[q[0] - 2, q[1] - 2], [q[0] + 2, q[1] + 2]]))
            xmax = abox.x1 - 0.025 * abox.width
            def clear(b):
                hit = ((dpx[:, 0] > b.x0 - rdot) & (dpx[:, 0] < b.x1 + rdot) & (dpx[:, 1] > b.y0 - rdot) & (dpx[:, 1] < b.y1 + rdot)).any()
                return (not hit) and b.x1 <= xmax and not any(b.overlaps(o) for o in blk)
            assert clear(bb), "slide label not clear at start: " + text
            dx = 0.0
            while clear(bb.translated(dx + 1.0, 0)): dx += 1.0
            stop = "right margin" if bb.translated(dx + 1.0, 0).x1 > xmax else "map dot / connector"
            bb = bb.translated(dx, 0); fx = fx + dx / pxu; t.set_position((fx, fy))
            # leader: shortest link from a box edge midpoint to the nearest point on the control segment
            pa_, pb_ = np.asarray(spec[2][0]), np.asarray(spec[2][1]); inv = ax.transData.inverted()
            def near(q):
                u = pb_ - pa_; tt = np.clip(np.dot(q - pa_, u) / np.dot(u, u), 0, 1); return pa_ + tt * u
            # candidates: points along the box top edge and right edge -> points along the segment; the leader may not
            # enter the zoom rectangle or any label box; keep the shortest
            rects = [r_.get_bbox() for r_ in ax.patches if isinstance(r_, mpl.patches.Rectangle) and not r_.get_fill()]
            tbox = [inv.transform(o_) for o_ in taken if o_.width > 20]   # placed labels / legend / inset (not dot markers)
            def blocked(q0, q1):
                for f_ in np.linspace(0.02, 1, 60):
                    q = q0 + f_ * (q1 - q0)
                    if any(r_.x0 <= q[0] <= r_.x1 and r_.y0 <= q[1] <= r_.y1 for r_ in rects): return True
                    if any(b_[0][0] <= q[0] <= b_[1][0] and b_[0][1] <= q[1] <= b_[1][1] for b_ in tbox): return True
                return False
            starts = [inv.transform((bb.x0 + f_ * bb.width, bb.y1)) for f_ in np.linspace(0.05, 0.95, 19)] + \
                     [inv.transform((bb.x1, bb.y0 + f_ * bb.height)) for f_ in np.linspace(0.2, 0.8, 4)]
            ends = [pa_ + f_ * (pb_ - pa_) for f_ in np.linspace(0, 1, 101)]
            cand = [(np.hypot(*(e_ - q_)), q_, e_) for q_ in starts for e_ in ends if not blocked(q_, e_)]
            assert cand, "no clear leader route"
            _, best, tgt = min(cand, key=lambda c_: c_[0])
            x, y = tgt; lead_from = best
            b_d = inv.transform(bb)
            if os.environ.get("FIG5_LAYOUT_CHECK"): print(f"LAYOUT_DI\ttext anchor (left, centre) x={fx:.4f} y={fy:.4f}\tbox x {b_d[0][0]:.4f}..{b_d[1][0]:.4f} "
                  f"y {b_d[0][1]:.4f}..{b_d[1][1]:.4f}\tslid {dx / pxu:.4f}\tstopped by {stop}\tright gap "
                  f"{(abox.x1 - bb.x1) / abox.width:.3f} of axis width\tleader ({best[0]:.4f}, {best[1]:.4f}) -> ({x:.4f}, {y:.4f})")
        else:
            bx = ax.transData.inverted().transform((bb.x1, (bb.y0 + bb.y1) / 2)); lead_from = (bx[0], fy)
        taken.append(bb)
        ax.annotate("", xy=(x, y), xytext=lead_from, arrowprops=dict(arrowstyle="-", lw=0.6, color="0.25", shrinkA=0, shrinkB=4), zorder=6)
    order = sorted([k for k in range(len(items)) if not isinstance(items[k][3], tuple)], key=lambda k: -items[k][0][1])
    for k in order:
        (x, y), text, _, side = items[k]; sgn = -1 if side == "left" else 1
        cands = []
        for r in (40, 60, 85, 115, 150, 190, 235):
            for ang in (0, 20, -20, 40, -40, 60, -60, 120, -120, 150, -150, 180):
                a = np.deg2rad(ang); dxp, dyp = sgn * r * np.cos(a), r * np.sin(a)
                cands.append((dxp, dyp))
        placed = None
        for dxp, dyp in cands:
            ha = "left" if dxp > 8 else ("right" if dxp < -8 else "center")
            t = ax.text(x + dxp / pxu, y + dyp / pxu, text, ha=ha, va="center", fontsize=fontsize, linespacing=1.15, zorder=7,
                        path_effects=[pe.withStroke(linewidth=3, foreground="white")])
            bb = t.get_window_extent(rr).expanded(1.03, 1.10)
            inside = abox.x0 + 2 < bb.x0 and bb.x1 < abox.x1 - 2 and abox.y0 + 2 < bb.y0 and bb.y1 < abox.y1 - 2
            if inside and not any(bb.overlaps(o) for o in taken):
                placed = (t, bb); break
            t.remove()
        assert placed is not None, text
        t, bb = placed; taken.append(bb)
        ax.annotate("", xy=(x, y), xytext=t.get_position(), arrowprops=dict(arrowstyle="-", lw=0.6, color="0.25", shrinkA=0, shrinkB=4), zorder=6)


def layout_panel_a(ax, Y, ix, clu_of, pairs, zoom, ins, legend, fontsize=9.5, report=False):
    """Panel a labels: one left-aligned column in the empty left margin, evenly spaced, ordered by the height of the
    point each leader attaches to, so that straight leaders cannot cross. Crohn disease - ulcerative colitis and
    psoriasis - psoriatic arthritis (all four in one cluster, markers within ~5 pt) share one label. The diabetes
    insipidus - T2D leader attaches to the end of the red line farthest from the zoom box. Returns the check table."""
    from scipy.spatial import ConvexHull
    from matplotlib.path import Path as MPath
    fig = ax.figure; fig.canvas.draw(); rr = fig.canvas.get_renderer(); tr = ax.transData; inv = tr.inverted()
    xl0, xl1 = ax.get_xlim(); yl0, yl1 = ax.get_ylim(); W, H = xl1 - xl0, yl1 - yl0
    zx0, zx1, zy0, zy1 = zoom
    P = {m: np.asarray(Y[ix[m]]) for p_ in pairs for m in (p_[1], p_[3])}
    cl = lambda m: int(clu_of.loc[m])
    def tail(a, b): return f"same cluster ({cl(a)})" if cl(a) == cl(b) else f"clusters {cl(a)} and {cl(b)}"
    byname = {p_[0]: p_ for p_ in pairs}
    nf, ts, di, cr, ps = (byname[k] for k in ("NF1 (mut)", "TSC", "diabetes insipidus", "Crohn disease", "psoriasis"))
    assert cl(cr[1]) == cl(cr[3]) == cl(ps[1]) == cl(ps[3])
    def far_end(p_):   # pair end farthest from the zoom box centre
        c = np.array([(zx0 + zx1) / 2, (zy0 + zy1) / 2]); return max((P[p_[1]], P[p_[3]]), key=lambda q: np.hypot(*(q - c)))
    def left_end(p_): return min((P[p_[1]], P[p_[3]]), key=lambda q: q[0])
    labs = [  # text, attach point, own pair keys (for the crossing check)
        (f"{nf[0]} \u2013 {nf[2]}\n{tail(nf[1], nf[3])}", left_end(nf), [nf]),
        (f"{ts[0]} \u2013 {ts[2]}\n{tail(ts[1], ts[3])}", left_end(ts), [ts]),
        (f"{di[0]} \u2013 {di[2]}\n{tail(di[1], di[3])}", far_end(di), [di]),
        (f"{cr[0]} \u2013\n{cr[2]};\n{ps[0]} \u2013\n{ps[2]};\nboth same cluster ({cl(cr[1])})", left_end(ps), [cr, ps]),
    ]
    labs.sort(key=lambda L: -L[1][1])
    # column: left edges on one line; top label level with its point; the DI label centred 0.09 above its point so
    # its near-horizontal leader clears the Crohn/psoriasis markers; the labels between and below are spaced with
    # equal vertical gaps (the gap between the top and DI labels sets the gap for all)
    xL = xl0 + 0.02 * W
    texts, boxes, leaders = [], [], []
    for txt, tgt, own in labs:
        texts.append(ax.text(xL, 0, txt, ha="left", va="center", fontsize=fontsize, linespacing=1.15, zorder=7,
                             path_effects=[pe.withStroke(linewidth=3, foreground="white")]))
    fig.canvas.draw()
    hts = [abs(inv.transform((0, t.get_window_extent(rr).y1))[1] - inv.transform((0, t.get_window_extent(rr).y0))[1]) for t in texts]
    kdi = [i for i, L in enumerate(labs) if L[0].startswith("diabetes insipidus")][0]
    c0, cdi = labs[0][1][1], labs[kdi][1][1] + 0.09
    span = (c0 - hts[0] / 2) - (cdi + hts[kdi] / 2) - sum(hts[1:kdi])
    gap = span / kdi; assert gap > 0.05 * H * 0 + 0.02 * H, gap
    yc = [c0]
    for i in range(1, len(labs)): yc.append(yc[-1] - hts[i - 1] / 2 - gap - hts[i] / 2)
    for t, y_ in zip(texts, yc): t.set_position((xL, y_))
    fig.canvas.draw()
    shr = 4 * fig.dpi / 72.0
    for t, (txt, tgt, own) in zip(texts, labs):
        bb = t.get_window_extent(rr); boxes.append(bb)
        s0 = np.array([bb.x1 + 3, (bb.y0 + bb.y1) / 2]); e0 = tr.transform(tgt)
        u = (e0 - s0) / np.hypot(*(e0 - s0)); e1 = e0 - u * shr
        ax.annotate("", xy=tuple(tgt), xytext=tuple(inv.transform(s0)),
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="0.25", shrinkA=0, shrinkB=4), zorder=6)
        leaders.append((s0, e1, own, txt))
    # ---------------- automatic check (display pixels) ----------------
    abox = ax.get_window_extent(rr); gx, gy = 0.02 * abox.width, 0.02 * abox.height
    zb = mpl.transforms.Bbox(tr.transform([[zx0, zy0], [zx1, zy1]]))
    dots = tr.transform(np.vstack([np.ma.getdata(c.get_offsets()) for c in ax.collections if len(c.get_offsets())]))
    hull = MPath(dots[ConvexHull(dots).vertices])
    pair_lines = [(tr.transform((l.get_xdata()[0], l.get_ydata()[0])), tr.transform((l.get_xdata()[1], l.get_ydata()[1])), l)
                  for l in ax.lines if len(l.get_xdata()) == 2]
    conns = []
    for a_ in ax.artists:
        if hasattr(a_, "_get_path_in_displaycoord"):
            v = a_._get_path_in_displaycoord().vertices; conns += [(v[i], v[i + 1]) for i in range(len(v) - 1)]
    def cross(p1, p2, q1, q2):
        d = lambda a, b, c: (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        return (d(p1, p2, q1) * d(p1, p2, q2) < 0) and (d(q1, q2, p1) * d(q1, q2, p2) < 0)
    def seg_box(p1, p2, b):
        return any(b.contains(*(p1 + f * (p2 - p1))) for f in np.linspace(0, 1, 200))
    def grown(b, fx, fy): return mpl.transforms.Bbox([[b.x0 - fx, b.y0 - fy], [b.x1 + fx, b.y1 + fy]])
    others = [legend.get_window_extent(rr), ins.get_tightbbox(rr), ax.title.get_window_extent(rr)]
    n = len(boxes); C = {}
    C["label - label (closer than 2%)"] = sum(grown(boxes[i], gx / 2, gy / 2).overlaps(grown(boxes[j], gx / 2, gy / 2)) for i in range(n) for j in range(i + 1, n))
    C["label - zoom box (closer than 2%)"] = sum(grown(b, gx, gy).overlaps(zb) for b in boxes)
    C["label - panel edge (closer than 2%)"] = sum(not (b.x0 >= abox.x0 + gx - 1 and b.x1 <= abox.x1 - gx and b.y0 >= abox.y0 + gy and b.y1 <= abox.y1 - gy) for b in boxes)
    C["label - legend / inset / title"] = sum(b.overlaps(o) for b in boxes for o in others)
    C["label - map dots"] = int(sum(((dots[:, 0] > b.x0) & (dots[:, 0] < b.x1) & (dots[:, 1] > b.y0) & (dots[:, 1] < b.y1)).sum() for b in boxes))
    C["label - dot hull"] = sum(hull.intersects_bbox(b) for b in boxes)
    C["leader - other label"] = sum(seg_box(s, e, boxes[j]) for i, (s, e, _, _) in enumerate(leaders) for j in range(n) if j != i)
    C["leader - leader crossing"] = sum(cross(leaders[i][0], leaders[i][1], leaders[j][0], leaders[j][1]) for i in range(n) for j in range(i + 1, n))
    C["leader - pair line crossing"] = sum(cross(s, e, q1, q2) for (s, e, _, _) in leaders for (q1, q2, _) in pair_lines)
    C["leader - inset connector crossing"] = sum(cross(s, e, q1, q2) for (s, e, _, _) in leaders for (q1, q2) in conns)
    C["leader - zoom box"] = sum(seg_box(s, e, zb) for (s, e, _, _) in leaders)
    dl = [L for L in leaders if L[3].startswith("diabetes insipidus")][0]
    mk = []
    for i, (s, e, own, _) in enumerate(leaders):
        ownm = {m for p_ in own for m in (p_[1], p_[3])}
        mk.append(sum(min(np.hypot(*((s + f * (e - s)) - tr.transform(P[m]))) for f in np.linspace(0, 1, 200)) < 8 for m in P if m not in ownm))
    C["leader - other pair marker (< 8 px)"] = sum(mk)
    C["DI leader within 3% of zoom box"] = int(seg_box(dl[0], dl[1], grown(zb, 0.03 * abox.width, 0.03 * abox.height)))
    if report:
        for k_, v_ in C.items(): print(f"CHECK_A\t{k_}\t{v_}")
        for t, b, (s, e, _, txt) in zip(texts, boxes, leaders):
            bd = inv.transform(b); sd, ed = inv.transform(s), inv.transform(e)
            print(f"LABEL_A\t{txt.splitlines()[0]}\tbox x {bd[0][0]:.3f}..{bd[1][0]:.3f} y {bd[0][1]:.3f}..{bd[1][1]:.3f}"
                  f"\tleader ({sd[0]:.3f}, {sd[1]:.3f}) -> ({ed[0]:.3f}, {ed[1]:.3f}) len {np.hypot(*(ed - sd)):.3f}\tgap {gap:.3f}")
    return C


def draw_panel_A_pairs(ax, pairs, singles, panel_letter, title_suffix,
                        marker_radius=0.03, min_visible_gap=0.10):
    # marker_radius: white-dot radius in data units at the current panel scale (s=24 marker,
    #   ~6.8in-wide axes over the padded data range -> ~0.024 data units radius).
    # min_visible_gap: minimum dash length we want showing *between* the two marker edges.
    #   Lines are only extended past the markers when the raw gap between the marker edges
    #   would otherwise be smaller than this (e.g. Wolman/CESD, whose markers nearly touch).
    present_count = 0
    x0, x1 = Y[:, 0].min(), Y[:, 0].max(); y0, y1 = Y[:, 1].min(), Y[:, 1].max()
    xc = (x0 + x1) / 2.0
    ax.scatter(Y[:, 0], Y[:, 1], c=point_colidx, cmap=cluster_cmap, vmin=0, vmax=n_clusters - 1,
               s=1.4, lw=0, rasterized=True, zorder=1)

    items = []
    for name_a, mid_a, name_b, mid_b, role in pairs:
        if mid_a not in ix or mid_b not in ix:
            print(f"panel {panel_letter} MISSING pair member: {name_a}/{name_b}")
            continue
        present_count += 1
        pa, pb = Y[ix[mid_a]], Y[ix[mid_b]]
        ca, cb = int(clu_of.loc[mid_a]), int(clu_of.loc[mid_b])
        color = GOOD if role == "pos" else BAD

        d = pb - pa
        raw_len = float(np.hypot(*d))
        gap_between_marker_edges = raw_len - 2 * marker_radius
        if gap_between_marker_edges < min_visible_gap and raw_len > 1e-9:
            u = d / raw_len
            ext = (min_visible_gap - gap_between_marker_edges) / 2.0
            la, lb = pa - u * ext, pb + u * ext
        else:
            la, lb = pa, pb
        line, = ax.plot([la[0], lb[0]], [la[1], lb[1]], ls="-", color=color, lw=2.0,
                         zorder=4, solid_capstyle="round", dash_capstyle="round")
        line.set_path_effects([pe.Stroke(linewidth=3.8, foreground="white"), pe.Normal()])
        for p in (pa, pb):
            ax.scatter(*p, s=24, facecolor="white", edgecolor="k", lw=0.8, zorder=5)

        tail = f"same cluster ({ca})" if ca == cb else f"clusters {ca} and {cb}"
        if role == "neg" and ca == cb:
            tail = f"same cluster ({ca})"
        name_line, _ = wrap_label_capped(f"{name_a} \u2013 {name_b}", start_width=24, max_lines=2)
        label_text = f"{name_line}\n{tail}"
        n_lines = label_text.count("\n") + 1
        mx, my = (pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0
        side = "left" if mx < xc else "right"
        if {mid_a, mid_b} == {"MONDO:0004782", "MONDO:0005148"}:   # DI-T2D: label in the empty lower-left margin
            side = ("slide", (x0 - 1.30, y0 + 0.10), (pa, pb))   # start lower-left, slid right in place_labels_local
        if {mid_a, mid_b} != {"MONDO:0016244", "MONDO:0100244"}:   # aHUS-PNH is labelled in the eculizumab inset
            items.append(((mx, my), label_text, n_lines, side))
    print(f"panel {panel_letter} present pairs: {present_count} / {len(pairs)}")
    for name_a, mid_a, name_b, mid_b, style in ECU_PAIRS:
        pa, pb = Y[ix[mid_a]], Y[ix[mid_b]]; ca, cb = int(clu_of.loc[mid_a]), int(clu_of.loc[mid_b])
        line, = ax.plot([pa[0], pb[0]], [pa[1], pb[1]], ls="-" if style == "solid" else (0, (1.2, 1.6)), color=GOOD,
                        lw=2.0, zorder=4, solid_capstyle="round", dash_capstyle="round")
        line.set_path_effects([pe.Stroke(linewidth=3.8, foreground="white"), pe.Normal()])
        for p in (pa, pb):
            ax.scatter(*p, s=24, facecolor="white", edgecolor="k", lw=0.8, zorder=5)
    # zoom inset on the four eculizumab-approved diseases (they sit too close together on the full map)
    E = {"aHUS": "MONDO:0016244", "PNH": "MONDO:0100244", "gMG": "MONDO:0009688", "NMOSD": "MONDO:0019100"}
    EP = np.array([Y[ix[m]] for m in E.values()]); pad = 0.06
    zx0, zx1 = EP[:, 0].min() - pad, EP[:, 0].max() + pad; zy0, zy1 = EP[:, 1].min() - pad, EP[:, 1].max() + pad
    ins = ax.inset_axes([1.04, 0.0, 0.36, 0.52])   # outside the map, right of panel a; ins.set_facecolor("white"); [sp.set_visible(True) or sp.set_color("0.4") for sp in ins.spines.values()]
    inw = (Y[:, 0] > zx0) & (Y[:, 0] < zx1) & (Y[:, 1] > zy0) & (Y[:, 1] < zy1)
    ins.scatter(Y[inw, 0], Y[inw, 1], s=3, c=cluster_cmap(point_colidx[inw]), alpha=0.5, lw=0, rasterized=True)
    for (na, ma, nb, mb) in [("aHUS", E["aHUS"], "PNH", E["PNH"], ), ("PNH", E["PNH"], "gMG", E["gMG"])]:
        ins.plot(*zip(Y[ix[ma]], Y[ix[mb]]), color=GOOD, lw=2.2, zorder=4)
    ins.plot(*zip(Y[ix[E["PNH"]]], Y[ix[E["NMOSD"]]]), color=GOOD, lw=2.2, ls=(0, (1.2, 1.6)), zorder=4)
    assert int(clu_of.loc[E["aHUS"]]) == int(clu_of.loc[E["PNH"]]) == int(clu_of.loc[E["gMG"]]) != int(clu_of.loc[E["NMOSD"]])
    for nm, m in E.items():
        p = Y[ix[m]]; yel = nm == "NMOSD"
        ins.scatter(*p, s=40, facecolor="#f2c400" if yel else "white", edgecolor="k", lw=0.9, zorder=5)
    ins_lab = [("PNH", "PNH", (-7, 9), "right"), ("aHUS", "aHUS", (7, 7), "left"), ("gMG", "gMG", (-9, 9), "right"), ("NMOSD", "NMOSD", (8, 0), "left")]
    for txt, key, off, ha in ins_lab:
        ins.annotate(f"{txt} ({int(clu_of.loc[E[key]])})", Y[ix[E[key]]], xytext=off, textcoords="offset points", ha=ha,
                     va="center", fontsize=9, zorder=6, path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
    ins.set_xlim(zx0, zx1); ins.set_ylim(zy0, zy1); ins.set_xticks([]); ins.set_yticks([])
    ins.set_title("eculizumab-approved diseases (cluster)", fontsize=9, pad=3,
                  path_effects=[pe.withStroke(linewidth=3, foreground="white")])
    from matplotlib.patches import Rectangle, ConnectionPatch
    ax.add_patch(Rectangle((zx0, zy0), zx1 - zx0, zy1 - zy0, fill=False, edgecolor="0.4", lw=0.8, zorder=5))
    for yy, iy in ((zy1, 1), (zy0, 0)):   # connectors from the right edge of the source box only (keeps clear of label leaders)
        ax.add_artist(ConnectionPatch(xyA=(zx1, yy), coordsA=ax.transData, xyB=(0, iy), coordsB=ins.transAxes,
                                      color="0.6", lw=0.7, zorder=4))
    for sp in ins.spines.values(): sp.set_visible(True); sp.set_edgecolor("0.4"); sp.set_linewidth(0.8)   # full frame (house style hides top/right)

    single_count = 0
    for name, mid in singles:
        if mid not in ix:
            print(f"panel {panel_letter} MISSING single: {name}")
            continue
        single_count += 1
        p = Y[ix[mid]]
        ax.scatter(*p, s=30, facecolor="#f2c400", edgecolor="k", lw=0.8, zorder=5)
        tagged = f"{name}\ncluster {int(clu_of.loc[mid])}"
        text, n_lines = tagged, 2
        side = "left" if p[0] < xc else "right"
        items.append(((p[0], p[1]), text, n_lines, side))
    print(f"panel {panel_letter} present singles: {single_count} / {len(singles)}")

    ax.set_aspect("equal"); ax.set_xlim(x0 - 1.35, x1 + 0.20); ax.set_ylim(y0 - 0.08, y1 + 0.02)  # same y range (scale) as panel b
    ax.xaxis.set_major_locator(MultipleLocator(0.6)); ax.yaxis.set_major_locator(MultipleLocator(0.3))
    ax.set_xlabel("Isomap 1", labelpad=3); ax.set_ylabel("Isomap 2", labelpad=2)
    ax.set_title(f"Diseases colored by cluster (K = {n_clusters}) \u2013 {title_suffix}", loc="left", fontsize=12)
    panel_label(ax, panel_letter)
    legA_handles = [Line2D([], [], color=GOOD, lw=2.0, label="positive control, same cluster"),
                    Line2D([], [], color=GOOD, lw=2.0, ls=(0, (1.2, 1.6)), label="positive control, neighboring cluster"),
                    Line2D([], [], color=BAD, lw=2.0, label="negative control"),
                    Line2D([], [], marker="o", color="none", markerfacecolor="#f2c400", markeredgecolor="k", markersize=7,
                           label="NMOSD (eculizumab-approved)")]
    leg = ax.legend(handles=legA_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False,
                    fontsize=9.5, borderpad=0.2, handlelength=1.8)   # outside the map, above the inset
    anchors = ([Y[ix[m]] for p_ in pairs for m in (p_[1], p_[3]) if m in ix] + [Y[ix[m]] for _, m in singles if m in ix]
               )
    assert not singles   # panel a has no single-point labels (NMOSD is labelled in the inset)
    layout_panel_a(ax, Y, ix, clu_of, pairs, (zx0, zx1, zy0, zy1), ins, leg, fontsize=9.5,
                   report=bool(os.environ.get("FIG5_LAYOUT_CHECK")))
    return legA_handles


def display_text(name):
    return name


def draw_panel_rarecommon(ax, targets, panel_letter, title_suffix):
    print(f"panel {panel_letter} pairs:", [(r_, int(clu_of.loc[a_])) for r_, a_, _, _ in targets])

    x0, x1 = Y[:, 0].min(), Y[:, 0].max(); y0, y1 = Y[:, 1].min(), Y[:, 1].max()
    is_rare = (cl["disease_class"].to_numpy() == "rare")
    # rare drawn first and semi-transparent; common drawn over it, opaque and slightly larger,
    # so the minority common class stays visible inside dense rare regions
    ax.scatter(Y[is_rare, 0], Y[is_rare, 1], c=RARE_COL, s=1.4, lw=0, alpha=0.35, rasterized=True, zorder=1, label="Rare")
    ax.scatter(Y[~is_rare, 0], Y[~is_rare, 1], c=COMMON_COL, s=2.0, lw=0, alpha=1.0, rasterized=True, zorder=2, label="Common")
    items, anchors = [], []
    xc = (x0 + x1) / 2.0
    for rname, ra, cname, ca_ in targets:
        pa, pb = Y[ix[ra]], Y[ix[ca_]]
        ln, = ax.plot([pa[0], pb[0]], [pa[1], pb[1]], color="0.15", lw=1.4, zorder=4, solid_capstyle="round")
        ln.set_path_effects([pe.Stroke(linewidth=3.2, foreground="white"), pe.Normal()])
        ax.scatter(*pa, s=30, facecolor="white", edgecolor=RARE_COL, lw=1.6, zorder=5)
        ax.scatter(*pb, s=30, facecolor="white", edgecolor=COMMON_COL, lw=1.6, zorder=5)
        anchors += [pa, pb]
        mx, my = (pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0
        side = "right" if pa[0] > 0.2 else "left"   # the right-hand example is labelled from the empty top-right corner
        rl = rname.split("\n")
        items.append(((pa, pb), tuple(rl) + (cname, f"cluster {int(clu_of.loc[ra])}"), (RARE_COL,) * len(rl) + (COMMON_COL, "0.35"), side))
    ax.set_aspect("equal"); ax.set_xlim(x0 - 1.35, x1 + 1.20); ax.set_ylim(y0 - 0.08, y1 + 0.02)  # same scale as panel a; right margin holds the legend
    ax.xaxis.set_major_locator(MultipleLocator(0.6)); ax.yaxis.set_major_locator(MultipleLocator(0.3))
    ax.set_xlabel("Isomap 1", labelpad=3)  # panel b: no y-label (shared with panel a)
    n_common, n_rare = int((~is_rare).sum()), int(is_rare.sum())
    ax.set_title(f"Common ({n_common:,}) and rare ({n_rare:,}) diseases \u2013 {title_suffix}", loc="left", fontsize=12)
    panel_label(ax, panel_letter)
    legB_handles = [Line2D([], [], marker="o", color="none", markerfacecolor=COMMON_COL, markersize=7, label="Common"),
                    Line2D([], [], marker="o", color="none", markerfacecolor=RARE_COL, markersize=7, label="Rare")]
    legB_handles += [Line2D([], [], marker="o", color="none", markerfacecolor="white", markeredgecolor=RARE_COL, mew=1.6,
                            markersize=7, label="rare example"),
                     Line2D([], [], marker="o", color="none", markerfacecolor="white", markeredgecolor=COMMON_COL, mew=1.6,
                            markersize=7, label="related common disease")]
    legB = ax.legend(handles=legB_handles, loc="lower right", bbox_to_anchor=(1.0, 0.0), frameon=False,
              fontsize=9.5, borderpad=0.3)   # inside the map, empty bottom-right corner
    import itertools
    fig = ax.figure; fig.canvas.draw(); rr = fig.canvas.get_renderer(); inv = ax.transData.inverted()
    XL = ax.get_xlim()[0] + 0.04; XR = ax.get_xlim()[1] - 0.04
    from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox
    # label lines: rare name in red, related common disease in blue, cluster in grey (colours carried per item)
    widths = []
    for _, lines, _, _ in items:   # widest line in data units (independent of the slot it is put in)
        w = 0
        for ln_ in lines:
            t = ax.text(XL, 0, ln_, ha="left", va="center", fontsize=9.5); bb = t.get_window_extent(rr); t.remove()
            w = max(w, inv.transform((bb.x1, 0))[0] - inv.transform((bb.x0, 0))[0])
        widths.append(w)
    SLOTS = [0.74, 0.34, -0.60, -1.00]   # left label centres: two above and two below the dense band left of the examples
    SLOTS_R = [0.98]                     # right label centre: empty top-right corner
    left = [k for k, it in enumerate(items) if it[3] == "left"]; right = [k for k, it in enumerate(items) if it[3] == "right"]
    assert len(left) <= len(SLOTS) and len(right) <= len(SLOTS_R), (len(left), len(right))
    segs = [it[0] for it in items]
    def _o(a, b, c): return np.sign((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]))
    def _x(p1, p2, p3, p4): return _o(p1, p2, p3) != _o(p1, p2, p4) and _o(p3, p4, p1) != _o(p3, p4, p2)
    def _d(p, a, b):
        a, b, p = map(np.asarray, (a, b, p)); u = b - a; t_ = np.clip(np.dot(p - a, u) / np.dot(u, u), 0, 1); return np.hypot(*(a + t_ * u - p))
    def _start(k, slot):
        return (XL + widths[k] + 0.03, slot) if items[k][3] == "left" else (XR - widths[k] - 0.03, slot)
    best = None
    for perm in itertools.permutations(range(len(SLOTS)), len(left)):
        slot_of = {k: SLOTS[perm[i]] for i, k in enumerate(left)}; slot_of.update({k: SLOTS_R[i] for i, k in enumerate(right)})
        for ends in itertools.product((0, 1), repeat=len(items)):
            L = [(_start(k, slot_of[k]), segs[k][ends[k]]) for k in range(len(items))]
            ok = all(not _x(*L[i], *L[j]) for i in range(len(L)) for j in range(i + 1, len(L)))
            ok = ok and all(not _x(*L[i], *segs[j]) for i in range(len(L)) for j in range(len(segs)) if j != i)
            ok = ok and all(_d(p, *L[i]) > 0.035 for i in range(len(L)) for j, sg in enumerate(segs) for p in sg if j != i)
            if ok:
                tot = sum(np.hypot(*(np.asarray(a) - np.asarray(b))) for a, b in L)
                if best is None or tot < best[0]: best = (tot, dict(slot_of), ends)
    assert best is not None, "no crossing-free label layout for panel b"
    _, slot_of, ends = best
    for k, (_, lines, cols, side) in enumerate(items):
        st = _start(k, slot_of[k])
        ax.annotate("", xy=segs[k][ends[k]], xytext=st, zorder=6,
                    arrowprops=dict(arrowstyle="-", color="0.3", lw=0.7, shrinkA=0, shrinkB=4))
        box = VPacker(children=[TextArea(t_, textprops=dict(color=c_, fontsize=9.5, path_effects=[pe.withStroke(linewidth=2.5, foreground="white")]))
                                for t_, c_ in zip(lines, cols)], align="left" if side == "left" else "right", pad=0, sep=1.5)
        ab = AnnotationBbox(box, (XL, slot_of[k]) if side == "left" else (XR, slot_of[k]), xycoords="data",
                            box_alignment=(0.0, 0.5) if side == "left" else (1.0, 0.5), frameon=False, pad=0, zorder=6)
        ax.add_artist(ab)
    perm = [None] * len(items)
    print("panel b label layout (slot, end):", [(slot_of[k], ends[k]) for k in range(len(items))])
    return legB_handles


# ---- Panels C/D: rarity-stratified transfer AUROC (All IxIDN-Matched pairs only; test-split dropped) + anchored ----
# Panel c comparators: extra score series drawn beside DisAtlas, with the same strata and the same bootstrap draws.
# Each entry: dict(key=<short id; output table OUTDIR/fig5_panel_c_<key>_auroc.tsv>, label=<legend text>,
#   scores=<TSV with blind_id and the score column, one row per benchmark pair, no missing values>, column=<score column,
#   higher = more likely co-investigated>, color, marker).
# TxGNN (published checkpoint), the knowledge-graph foundation model, with the scores and protocol of the cross-system
# figure (fig5c_export_txgnn_scores.py; its 213 unmappable pairs ranked last). DisAtlas keeps the orange of this figure;
# TxGNN is dark brown #5c3317 (luminance contrast 3.7:1 against the orange).
PANEL_C_COMPARATORS = [
    dict(key="txgnn", label="TxGNN", scores=os.path.join(DV, "fig5_clustering", "fig5c_txgnn_pair_scores.tsv"), column="score",
         color="#5c3317", marker="^"),
]
PANEL_C_DISATLAS_COLOR = "#d9822b"   # DisAtlas colour of this figure
assert len({c_["key"] for c_ in PANEL_C_COMPARATORS}) == len(PANEL_C_COMPARATORS), "duplicate comparator key"

# ---- Panel D data: DisAtlas score p(i,j) by stratum on the 4,364 IxIDN-Matched pairs (y=1 true pair, y=0 non-pair).
S_SCORES = pd.read_csv(os.path.join(DV, "fig5_clustering", "ixidn_matched_pair_scores.tsv"), sep="\t")
assert len(S_SCORES) == 4364 and S_SCORES.a.isin(ix).all() and S_SCORES.b.isin(ix).all()
# the figure reads p back from the released distance matrix (d = 1 - p, float32): p = 1 - d in float32
S_SCORES["p"] = (np.float32(1) - (np.float32(1) - S_SCORES.p.to_numpy(np.float32))).astype(float)
from sklearn.metrics import roc_auc_score
# 95% CIs: disease-cluster bootstrap (each disease with all its pairs is one resampling cluster, pair weight w_a*w_b),
# one set of draws shared by DisAtlas, every comparator series and every stratum, and by the SI TxGNN comparison and the
# SI rarity strata: disease order = pd.unique over (a, b) of the cross-system per-run table
# (data/derived/fig6_comparison/per_run_scores.tsv), 10 shards of B/10 draws, shard s seeded np.random.default_rng([4707, s]).
_B=int(os.environ.get("FIG5_B","2000")); assert _B%10==0, "FIG5_B must be a multiple of the 10 shards"
for _cmp in PANEL_C_COMPARATORS:
    _W=pd.read_csv(_cmp["scores"],sep="\t"); _col=_cmp["column"]
    _m=S_SCORES[["blind_id"]].merge(_W[["blind_id",_col]],on="blind_id",how="left",validate="one_to_one")
    assert len(_m)==len(S_SCORES) and _m[_col].notna().all(), ("comparator must score every benchmark pair", _cmp["key"])
    S_SCORES["cmp_"+_cmp["key"]]=_m[_col].to_numpy(float)
_SERIES=["p"]+["cmp_"+_cmp["key"] for _cmp in PANEL_C_COMPARATORS]
_PR=pd.read_csv(os.path.join(DV, "fig6_comparison", "per_run_scores.tsv"),sep="\t",usecols=["a","b"]); _nodes=pd.unique(_PR[["a","b"]].values.ravel())
assert set(_nodes)==set(pd.unique(S_SCORES[["a","b"]].values.ravel())), "cluster set differs from the SI bootstrap"
assert set(map(frozenset,_PR[["a","b"]].values))==set(map(frozenset,S_SCORES[["a","b"]].values)), "pair set differs"
_WB=[]
for _s in range(10):
    _rg=np.random.default_rng([4707,_s]); _WB+=[np.bincount(_rg.integers(0,len(_nodes),len(_nodes)),minlength=len(_nodes)) for _ in range(_B//10)]
_nix={d_:i for i,d_ in enumerate(_nodes)}
_ai=S_SCORES.a.map(_nix).to_numpy(); _bi=S_SCORES.b.map(_nix).to_numpy(); print("panel c bootstrap B",len(_WB),"clusters",len(_nodes))
_res={c_:[] for c_ in _SERIES}; _diff=[]
for _st in ["rare-rare","rare-common","common-common","all"]:
    _k=np.ones(len(S_SCORES),bool) if _st=="all" else (S_SCORES.rare_common_pair==_st).to_numpy()
    y_=S_SCORES.y.to_numpy()[_k]; sc_={c_:S_SCORES[c_].to_numpy()[_k] for c_ in _SERIES}; bs={c_:[] for c_ in _SERIES}
    for w_ in _WB:
        ww_=(w_[_ai]*w_[_bi])[_k].astype(float); q_=ww_>0
        if len(set(y_[q_]))<2: continue
        for c_ in bs: bs[c_].append(roc_auc_score(y_[q_],sc_[c_][q_],sample_weight=ww_[q_]))
    for c_ in _SERIES[1:]:
        _d=np.asarray(bs["p"])-np.asarray(bs[c_])
        _diff.append(dict(stratum=_st,comparator=c_[4:],disatlas=roc_auc_score(y_,sc_["p"]),comparator_auroc=roc_auc_score(y_,sc_[c_]),
                          diff=roc_auc_score(y_,sc_["p"])-roc_auc_score(y_,sc_[c_]),ci_lo=np.percentile(_d,2.5),ci_hi=np.percentile(_d,97.5),
                          n_boot=len(_d)))
    for c_ in _SERIES:
        _res[c_].append(dict(stratum=_st,auroc=roc_auc_score(y_,sc_[c_]),ci_lo=np.percentile(bs[c_],2.5),ci_hi=np.percentile(bs[c_],97.5),
                        n_pos=int((y_==1).sum()),n_neg=int((y_==0).sum()),n_boot=len(bs[c_])))
AUC32=pd.DataFrame(_res["p"]).set_index("stratum"); AUC32.to_csv(f"{OUT}/fig5_panel_c_auroc.tsv",sep="\t"); print(AUC32)
CMP_AUC={}
PAIRED_DIFF=pd.DataFrame(_diff)
if PANEL_C_COMPARATORS:
    PAIRED_DIFF.to_csv(f"{OUT}/fig5_panel_c_paired_differences.tsv",sep="\t",index=False); print(PAIRED_DIFF)
for _cmp in PANEL_C_COMPARATORS:
    CMP_AUC[_cmp["key"]]=pd.DataFrame(_res["cmp_"+_cmp["key"]]).set_index("stratum")
    CMP_AUC[_cmp["key"]].to_csv(f"{OUT}/fig5_panel_c_{_cmp['key']}_auroc.tsv",sep="\t"); print(_cmp["key"]); print(CMP_AUC[_cmp["key"]])
assert {"a", "b", "p", "y", "rare_common_pair"}.issubset(S_SCORES.columns)


def _logit(p, eps=1e-12):
    p = np.clip(np.asarray(p, dtype=np.float64), eps, 1 - eps)
    return np.log(p / (1 - p))


SCORE_DIST = {}
for _st in ["rare-rare", "rare-common", "common-common", "all"]:
    _d = S_SCORES if _st == "all" else S_SCORES[S_SCORES.rare_common_pair == _st]
    SCORE_DIST[_st] = {"pos": _logit(_d.loc[_d.y == 1, "p"].to_numpy()),
                        "neg": _logit(_d.loc[_d.y == 0, "p"].to_numpy())}

# Control pairs (panel A) are a separate set, not part of the IxIDN-Matched benchmark; their p(i,j) comes from the
# released pair-score matrix and their stratum is the pair's own rarity classification.
CONTROL_PAIR_SCORES = []
for _nameA, _idA, _nameB, _idB, _role in PANEL_A_PAIRS:
    _p_ctrl = float(PAIR_P[(_idA, _idB)])
    _rareA = bool(cl.loc[_idA, "is_rare"])
    _rareB = bool(cl.loc[_idB, "is_rare"])
    _stratum = "rare-rare" if (_rareA and _rareB) else ("common-common" if not (_rareA or _rareB) else "rare-common")
    CONTROL_PAIR_SCORES.append(dict(name=f"{_nameA} / {_nameB}", p=_p_ctrl, logit=_logit([_p_ctrl])[0],
                                     role=_role, stratum=_stratum))
for _d in CONTROL_PAIR_SCORES:
    print("control", _d["name"], _d["role"], _d["stratum"], round(_d["p"], 4))
    assert not (_d["role"] == "pos" and _d["p"] < 0.5), ("positive control below 0.5 -- report before rendering", _d)
ST = [("rare-rare", "Rare-\nrare"), ("rare-common", "Rare-\ncommon"), ("common-common", "Common-\ncommon"), ("all", "All\npairs")]


# Stratum pair counts (benchmark design; common-common has only 26 negatives -> hollow marker, low-precision estimate)
N_EXP = {"rare-rare": (1284, 1168), "rare-common": (253, 988), "common-common": (645, 26), "all": (2182, 2182)}
HOLLOW_C = {"common-common"}


def draw_panel_c_allpairs(ax, panel_letter):
    assert [st for st, _ in ST] == list(N_EXP)
    DIS_COL = PANEL_C_DISATLAS_COLOR
    # series: DisAtlas first, then the configured comparators; offsets spread symmetrically around each stratum
    series = [dict(tab=AUC32, label="DisAtlas", color=DIS_COL, marker="o", ms=6, lw=1.4, mew=1.4)] + \
             [dict(tab=CMP_AUC[c_["key"]], label=c_["label"], color=c_.get("color", "#6b7f99"), marker=c_.get("marker", "s"),
                   ms=5.5, lw=1.2, mew=1.3) for c_ in PANEL_C_COMPARATORS]
    offs = (np.arange(len(series)) - (len(series) - 1) / 2) * 0.18
    # lower axis limit: 0.30 unless a CI reaches below it (then the next 0.05 step under the lowest CI)
    ymin = min(0.30, np.floor((min(s_["tab"].ci_lo.min() for s_ in series) - 0.02) * 20) / 20)
    for i, (st, lab) in enumerate(ST):
        hol = st in HOLLOW_C
        for s_, dx in zip(series, offs):
            q = s_["tab"].loc[st]
            assert (int(q.n_pos), int(q.n_neg)) == N_EXP[st], (st, s_["label"])
            ax.errorbar(i + dx, q.auroc, yerr=[[q.auroc - q.ci_lo], [q.ci_hi - q.auroc]], fmt=s_["marker"], color=s_["color"],
                        ms=s_["ms"], capsize=3, lw=s_["lw"], zorder=3, label=s_["label"] if i == 0 else None,
                        mfc="white" if hol else s_["color"], mec=s_["color"], mew=s_["mew"])
        q = AUC32.loc[st]
        ax.text(i, ymin + 0.03, f"{int(q.n_pos)}/{int(q.n_neg)}", ha="center", fontsize=9.5, color="0.35", zorder=5,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0))
        if hol:
            top = max(s_["tab"].loc[st].ci_hi for s_ in series)
            ax.text(i, top + 0.025, f"hollow = few\nnegatives (n = {int(q.n_neg)})", ha="center", va="bottom",
                    fontsize=9, color="0.35", linespacing=1.1)
    ax.axhline(0.5, color="0.8", ls="--", lw=0.9)
    ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l in ST])
    ax.set_ylim(ymin, 1.08); ax.set_xlim(-0.5, 3.5)
    ax.set_ylabel("AUROC (95% CI)")
    ax.set_title("All benchmark pairs, by pair rarity", loc="left", fontsize=12)
    panel_label(ax, panel_letter)
    ax.legend(loc="upper left", frameon=False, fontsize=10.5, bbox_to_anchor=(0.0, 1.02), ncol=len(series), handletextpad=0.3, columnspacing=1.0)


def draw_panel_d_scoredist(ax, panel_letter):
    """UPNA-DDAL score p(i,j) distributions (logit axis), true-pair vs non-pair, by rarity stratum.
    Overlaid: panel A's 5 validated control pairs as colored ticks on their own rarity stratum."""
    ST_D = [("rare-rare", "Rare-\nrare"), ("rare-common", "Rare-\ncommon"),
            ("common-common", "Common-\ncommon"), ("all", "All\npairs")]
    assert [x for x, _ in ST_D] == [x for x, _ in ST] == list(N_EXP)
    OFF, WID = 0.18, 0.30
    parts_pos, parts_neg = [], []
    for i, (st, lab) in enumerate(ST_D):
        for sign, off, color in [("pos", -OFF, GOOD), ("neg", OFF, BAD)]:
            vals = SCORE_DIST[st][sign]
            vp = ax.violinplot([vals], positions=[i + off], widths=WID, showmeans=False,
                                showmedians=False, showextrema=False)
            ax.scatter([i + off], [np.median(vals)], s=26, marker="o", facecolor="white", edgecolor="0.15", lw=1.0, zorder=6)
            for b in vp["bodies"]:
                b.set_facecolor(color); b.set_edgecolor(color); b.set_alpha(0.42); b.set_zorder(2)
            (parts_pos if sign == "pos" else parts_neg).append(vp["bodies"][0])
        npos, nneg = len(SCORE_DIST[st]["pos"]), len(SCORE_DIST[st]["neg"])
        assert (npos, nneg) == N_EXP[st], (st, npos, nneg)   # same strata order and counts as panel c
        ax.text(i, -11.1, f"{npos}+ / {nneg}\u2212", ha="center", fontsize=9.5, color="0.35")
    # control pairs: tick at own p(i,j), on the x-offset matching its pos/neg role, within its own
    # rarity stratum (verified against disease_meta.is_rare, not inferred from name/role)
    stx = {st: i for i, (st, _) in enumerate(ST_D)}
    for d in CONTROL_PAIR_SCORES:
        i = stx[d["stratum"]]
        off = -OFF if d["role"] == "pos" else OFF
        color = GOOD if d["role"] == "pos" else BAD
        ax.plot([i + off - WID * 0.42, i + off + WID * 0.42], [d["logit"], d["logit"]],
                 color=color, lw=2.6, solid_capstyle="butt", zorder=5,
                 path_effects=[pe.Stroke(linewidth=4.0, foreground="white"), pe.Normal()])
    ax.axhline(0.0, color="0.8", ls="--", lw=0.9, zorder=1)
    ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l in ST_D])
    ax.set_xlim(-0.55, 3.55); ax.set_ylim(-11.6, 7.6)
    ptick = [0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999]
    ax.set_yticks(_logit(ptick)); ax.set_yticklabels([str(x) for x in ptick])
    ax.set_ylabel(r"DisAtlas score  $p_{i,j}$  (logit scale)")
    ax.set_title("Score distribution by pair rarity", loc="left", fontsize=12)
    panel_label(ax, panel_letter)
    parts_pos[0].set_label("True pair ($y=1$)"); parts_neg[0].set_label("Non-pair ($y=0$)")
    leg_handles = [parts_pos[0], parts_neg[0],
                   Line2D([], [], color=GOOD, lw=2.6, solid_capstyle="butt", label="positive control pair"),
                   Line2D([], [], color=BAD, lw=2.6, solid_capstyle="butt", label="negative control pair"),
                   Line2D([], [], marker="o", ls="", mfc="white", mec="0.15", ms=5.5, label="median")]
    ax.legend(handles=leg_handles, loc="upper center", ncol=3, fontsize=9.5,
              handlelength=1.6, columnspacing=1.1, labelspacing=0.35, handletextpad=0.5,
              framealpha=0.9, edgecolor="none", borderpad=0.6)

fig = plt.figure(figsize=(17.7, 9.6))   # top row: map | inset + legend | map (legend inside); same map scale in a and b
axA = fig.add_axes([0.0158, 0.48803, 0.3406, 0.47587])
axB = fig.add_axes([0.5650, 0.48803, 0.4236, 0.47587])
axC = fig.add_axes([0.02, 0.1052, 0.43165, 0.3497])
axD = fig.add_axes([0.54165, 0.1052, 0.43165, 0.3497])

draw_panel_A_pairs(axA, PANEL_A_PAIRS, PANEL_A_SINGLES, "a", "control pairs")
draw_panel_rarecommon(axB, PANEL_B_PAIRS, "b", "rare examples")
draw_panel_c_allpairs(axC, "c")
draw_panel_d_scoredist(axD, "d")

# write the figure straight to the manuscript FIGURES folder (no stray file in the working directory)
for _out in (f"{OUT}/disease_clustering_map.png", f"{OUT}/disease_clustering_map.pdf"):
    fig.savefig(_out, dpi=int(os.environ.get("FIG_DPI", "280")), bbox_inches="tight")  # 280 dpi keeps the PNG under 1.9 MB
    print("saved", _out)

if os.environ.get("FIG5_AB"):   # optional review render of panels a and b alone (not used by the build)
  fig_ab = plt.figure(figsize=(17.7, 5.2))
  a2 = fig_ab.add_axes([0.0158, 0.05, 0.3406, 0.88]); b2 = fig_ab.add_axes([0.5650, 0.05, 0.4236, 0.88])
  draw_panel_A_pairs(a2, PANEL_A_PAIRS, PANEL_A_SINGLES, "a", "control pairs")
  draw_panel_rarecommon(b2, PANEL_B_PAIRS, "b", "rare examples")
  for _out in (f"{OUT}/fig5_disease_clustering_ab.png", f"{OUT}/fig5_disease_clustering_ab.pdf"):
    fig_ab.savefig(_out, dpi=int(os.environ.get("FIG_DPI", "280")), bbox_inches="tight"); print("saved", _out)

from disatlas.parity import write_values
_V = {"n_diseases": n, "n_clusters": 55}
for _st in AUC32.index:
    for _c_ in ("auroc", "ci_lo", "ci_hi", "n_pos", "n_neg"):
        _V[f"disatlas|{_st}|{_c_}"] = AUC32.loc[_st, _c_]
        for _k_, _t_ in CMP_AUC.items(): _V[f"{_k_}|{_st}|{_c_}"] = _t_.loc[_st, _c_]
for _r_ in PAIRED_DIFF.itertuples():
    for _c_ in ("diff", "ci_lo", "ci_hi"): _V[f"diff_{_r_.comparator}|{_r_.stratum}|{_c_}"] = getattr(_r_, _c_)
for _d in CONTROL_PAIR_SCORES: _V[f"control|{_d['name']}"] = round(_d["p"], 6)
write_values(f"{OUT}/disease_clustering_map_values.tsv", _V)
