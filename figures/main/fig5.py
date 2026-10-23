"""Figure 5 (Results, "DisAtlas organises diseases by shared mechanism"): disease map and rarity-stratified transfer.

  a  five validated control pairs on the K = 50 UPNA-DDAL cluster map (weighted-MDS display coordinates)
  b  common / rare disease map with named clinical examples
  c  IxIDN-Matched AUROC by pair-rarity stratum (DisAtlas with 95% node-bootstrap CI; HPO and GO Resnik baselines)
  d  UPNA-DDAL score p(i,j) distributions, true pair vs non-pair, by stratum; panel a control pairs as ticks

Inputs (repo-relative):
  figures/main/data/fig5_display_map.tsv.gz         disease_id, cluster, x, y (13,019 diseases)
  figures/main/data/fig5_control_pairs.tsv          the five control pairs and their UPNA-DDAL score p
  data/08S_eval/ixidn4364_scores.tsv                IxIDN-Matched pairs: y, rarity stratum, DisAtlas probability p
  data/00_rebuilt_ontology/ixidn_matched_hpo_scores.csv  HPO Resnik-BMA per pair (rebuilt by code/00_ontology/)
  figures/main/data/ixidn_matched_go_resnik.tsv.gz  GO Resnik gene-BMA per pair
  data/00_raw/llm_annotations/disease_labels.parquet is_rare per disease (stored labels; no model is called)

    python figures/main/fig5.py      -> figures/main/fig5_disease_map.{png,pdf}
"""
import os, sys, textwrap
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from disatlas_style import apply_style, panel_label, save, REPO_ROOT, FIG_DATA
from fig_metrics import auroc, node_bootstrap_auroc

cl = pd.read_csv(os.path.join(FIG_DATA, "fig5_display_map.tsv.gz"), sep="\t").rename(columns={"x": "mds1", "y": "mds2"})
assert len(cl) == 13019 and cl["mds1"].notna().all() and cl["mds2"].notna().all()
dm = pd.read_parquet(os.path.join(REPO_ROOT, "data/00_raw/llm_annotations/disease_labels.parquet")).set_index("disease_id")
cl = cl.set_index("disease_id").join(dm[["disease_class", "is_rare"]], how="left")
assert cl["disease_class"].notna().all()
assert (int((cl["is_rare"] == 1).sum()), int((cl["is_rare"] == 0).sum())) == (10488, 2531)  # Results and Figure 5 caption
_net = pd.read_parquet(os.path.join(REPO_ROOT, "data/03_network/mechdist_network_nodes.parquet.gzip")).disease_id.map(dm.is_rare)
assert (len(_net), int((_net == 1).sum()), int((_net == 0).sum())) == (13021, 10489, 2532)  # Methods, network coverage
ids = cl.index.to_numpy(); Y = cl[["mds1", "mds2"]].to_numpy(); ix = {d: i for i, d in enumerate(ids)}
clu_of = cl["cluster"]

# ---- Panel A: 5 validated control pairs (K = 50 UPNA-DDAL clusters) ----
GOOD, BAD = "#00b140", "#b42318"
PANEL_A_PAIRS = [
    ("Wolman (LAL-D)", "MONDO:0019148", "CESD", "MONDO:0019149", "pos"),
    ("aHUS", "MONDO:0035290", "PNH", "MONDO:0100244", "pos"),
    ("Crohn disease", "MONDO:0005011", "ulcerative colitis", "MONDO:0005101", "pos"),
    ("diabetes insipidus", "MONDO:0004782", "T2D", "MONDO:0005148", "neg"),
    ("NF1 (mut)", "MONDO:0018208", "NF2", "MONDO:0007039", "neg"),
]
PANEL_A_SINGLES = []

PANEL_B_TARGETS = [
    ("Infantile dystonia-parkinsonism", "MONDO:0013150", "left"),
    ("DYT17", "MONDO:0012895", "left"),
    ("Dystonia 16", "MONDO:0012789", "right"),
    ("ET5", "MONDO:0014756", "right"),
    ("PARK5", "DIS_B12E133F200", "right"),
    ("HMPS", "MONDO:0011023", "right"),
    ("LFS", "MONDO:0018875", "left"),
    ("Melanoma", "MONDO:0005105", "right"),
    ("Prostate carcinoma", "MONDO:0005159", "left"),
    ("Dystonia-aphonia syndrome", "MONDO:0018476", None),
]

import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import patheffects as pe
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import ListedColormap

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
                               push_to_edge_frac=0.0):
    y_mid = (y0 + y1) / 2.0
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
            if region == "bottom" and lo_local < y0 + edge_gap:
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


def draw_panel_A_pairs(ax, pairs, singles, panel_letter, title_suffix,
                        marker_radius=0.024, min_visible_gap=0.045):
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
        line, = ax.plot([la[0], lb[0]], [la[1], lb[1]], ls=(0, (1, 1.3)), color=color, lw=2.6,
                         zorder=4, solid_capstyle="round", dash_capstyle="round")
        line.set_path_effects([pe.Stroke(linewidth=4.4, foreground="white"), pe.Normal()])
        for p in (pa, pb):
            ax.scatter(*p, s=24, facecolor="white", edgecolor="k", lw=0.8, zorder=5)

        tail = f"shared cluster (c{ca})" if ca == cb else f"separate clusters (c{ca}, c{cb})"
        if role == "neg" and ca == cb:
            tail = f"shares cluster (c{ca})"
        name_line, _ = wrap_label_capped(f"{name_a} - {name_b}", start_width=24, max_lines=2)
        label_text = f"{name_line}\n{tail}"
        n_lines = label_text.count("\n") + 1
        mx, my = (pa[0] + pb[0]) / 2.0, (pa[1] + pb[1]) / 2.0
        side = "left" if mx < xc else "right"
        items.append(((mx, my), label_text, n_lines, side))
    print(f"panel {panel_letter} present pairs: {present_count} / {len(pairs)}")

    single_count = 0
    for name, mid in singles:
        if mid not in ix:
            print(f"panel {panel_letter} MISSING single: {name}")
            continue
        single_count += 1
        p = Y[ix[mid]]
        ax.scatter(*p, s=24, facecolor="white", edgecolor="k", lw=0.8, zorder=5)
        tagged = f"{name} (c{int(clu_of.loc[mid])})"
        text, n_lines = wrap_label_capped(tagged, start_width=20, max_lines=2)
        side = "left" if p[0] < xc else "right"
        items.append(((p[0], p[1]), text, n_lines, side))
    print(f"panel {panel_letter} present singles: {single_count} / {len(singles)}")

    place_leader_labels_sided(ax, items, x0, x1, y0, y1, fontsize=10, line_height_frac=0.068, push_to_edge_frac=0.6)
    ax.set_aspect("equal"); ax.set_xlim(x0 - 1.35, x1 + 0.95); ax.set_ylim(y0 - 0.08, y1 + 0.02)
    ax.xaxis.set_major_locator(MultipleLocator(0.3)); ax.yaxis.set_major_locator(MultipleLocator(0.3))
    ax.set_xlabel("Weighted MDS 1", labelpad=3); ax.set_ylabel("Weighted MDS 2", labelpad=2)
    ax.set_title(f"Diseases colored by cluster assignment (k={n_clusters}) -- {title_suffix}", loc="left", fontsize=12)
    panel_label(ax, panel_letter)
    legA_handles = [Line2D([], [], color=GOOD, lw=1.6, ls=(0, (1.5, 1.8)), label="positive\ncontrol pair"),
                    Line2D([], [], color=BAD, lw=1.6, ls=(0, (1.5, 1.8)), label="negative\ncontrol pair")]
    ax.legend(handles=legA_handles, loc="lower right", bbox_to_anchor=(1.0, 0.0), frameon=True,
              framealpha=0.9, edgecolor="none", fontsize=11, borderpad=0.6)
    return legA_handles


def display_text(name):
    return name


def draw_panel_rarecommon(ax, targets, panel_letter, title_suffix):
    present, missing = [], []
    for name, mid, side in targets:
        (present if mid in ix else missing).append((name, mid, side) if mid in ix else (name, mid))
    print(f"panel {panel_letter} present:", present)
    print(f"panel {panel_letter} missing:", missing)

    x0, x1 = Y[:, 0].min(), Y[:, 0].max(); y0, y1 = Y[:, 1].min(), Y[:, 1].max()
    is_rare = (cl["disease_class"].to_numpy() == "rare")
    # rare drawn first and semi-transparent; common drawn over it, opaque and slightly larger,
    # so the minority common class stays visible inside dense rare regions (author request)
    ax.scatter(Y[is_rare, 0], Y[is_rare, 1], c=RARE_COL, s=1.4, lw=0, alpha=0.35, rasterized=True, zorder=1, label="Rare")
    ax.scatter(Y[~is_rare, 0], Y[~is_rare, 1], c=COMMON_COL, s=2.0, lw=0, alpha=1.0, rasterized=True, zorder=2, label="Common")
    items = []
    for name, mid, side in present:
        p = Y[ix[mid]]
        ax.scatter(*p, s=24, facecolor="white", edgecolor="k", lw=0.8, zorder=5)
        tagged = f"{display_text(name)} (c{int(clu_of.loc[mid])})"
        text, n_lines = wrap_label_capped(tagged)
        items.append(((p[0], p[1]), text, n_lines, side))
    place_leader_labels_sided(ax, items, x0, x1, y0, y1, fontsize=10)
    ax.set_aspect("equal"); ax.set_xlim(x0 - 1.35, x1 + 0.95); ax.set_ylim(y0 - 0.08, y1 + 0.02)  # unified w/ panel A: same scale/ticks, no distortion
    ax.xaxis.set_major_locator(MultipleLocator(0.3)); ax.yaxis.set_major_locator(MultipleLocator(0.3))
    ax.set_xlabel("Weighted MDS 1", labelpad=3)  # panel b: no y-label (shared with panel a)
    n_common, n_rare = int((~is_rare).sum()), int(is_rare.sum())
    ax.set_title(f"Common ({n_common:,}) and rare ({n_rare:,}) diseases -- {title_suffix}", loc="left", fontsize=12)
    panel_label(ax, panel_letter)
    legB_handles = [Line2D([], [], marker="o", color="none", markerfacecolor=COMMON_COL, markersize=7, label="Common"),
                    Line2D([], [], marker="o", color="none", markerfacecolor=RARE_COL, markersize=7, label="Rare")]
    ax.legend(handles=legB_handles, loc="center right", bbox_to_anchor=(1.0, 0.5), frameon=True,
              framealpha=0.9, edgecolor="none", fontsize=11, borderpad=0.6)
    return legB_handles


# ---- Panels C/D data: IxIDN-Matched pairs; DisAtlas column = UPNA-DDAL probability p(i,j); y=1 true pair ----
S_SCORES = pd.read_csv(os.path.join(REPO_ROOT, "data/08S_eval/ixidn4364_scores.tsv"), sep="\t",
                       usecols=["blind_id", "a", "b", "y", "p", "rare_common_pair"])
S_SCORES = S_SCORES.merge(pd.read_csv(os.path.join(REPO_ROOT, "data/00_rebuilt_ontology/ixidn_matched_hpo_scores.csv"))[
    ["blind_id", "hpo_resnik_bma"]], on="blind_id", how="left")
S_SCORES = S_SCORES.merge(pd.read_csv(os.path.join(FIG_DATA, "ixidn_matched_go_resnik.tsv.gz"), sep="\t"), on="blind_id", how="left")
assert len(S_SCORES) == 4364 and int(S_SCORES.y.sum()) == 2182
# rarity stratum re-derived from the stored disease labels and checked against the shipped column
_ra = S_SCORES.a.map(dm.is_rare).astype(int); _rb = S_SCORES.b.map(dm.is_rare).astype(int)
_rc = np.where((_ra == 1) & (_rb == 1), "rare-rare", np.where((_ra + _rb) == 1, "rare-common", "common-common"))
assert (_rc == S_SCORES.rare_common_pair).all()
_pen = lambda s: s.fillna(s.min() - 1)   # unscorable pairs ranked last (over all 4,364 pairs)
S_SCORES["HPO Resnik-BMA"] = _pen(S_SCORES.hpo_resnik_bma).astype(float)
S_SCORES["GO Resnik gene-BMA"] = _pen(S_SCORES.go_resnik_genebma).astype(float)
_rows = []
for _st in ("rare-rare", "rare-common", "common-common", "all"):
    _D = S_SCORES if _st == "all" else S_SCORES[S_SCORES.rare_common_pair == _st]
    _bs = node_bootstrap_auroc(_D.a.values, _D.b.values, _D.y.values, _D.p.values, B=1000, seed=4707)
    _rows.append({"stratum": _st, "auroc": auroc(_D.y.values, _D.p.values),
                  "ci_lo": np.percentile(_bs, 2.5), "ci_hi": np.percentile(_bs, 97.5),
                  "HPO Resnik-BMA | auroc": auroc(_D.y.values, _D["HPO Resnik-BMA"].values),
                  "GO Resnik gene-BMA | auroc": auroc(_D.y.values, _D["GO Resnik gene-BMA"].values),
                  "n_pos": int(_D.y.sum()), "n_neg": int((1 - _D.y).sum())})
T = pd.DataFrame(_rows); T["scope"] = "All 4364 pairs"
T.to_csv(os.path.join(HERE, "fig5c_stratified_auroc.tsv"), sep="\t", index=False, float_format="%.6f")


def _logit(p, eps=1e-12):
    p = np.clip(np.asarray(p, dtype=np.float64), eps, 1 - eps)
    return np.log(p / (1 - p))


SCORE_DIST = {}
for _st in ["rare-rare", "rare-common", "common-common", "all"]:
    _d = S_SCORES if _st == "all" else S_SCORES[S_SCORES.rare_common_pair == _st]
    SCORE_DIST[_st] = {"pos": _logit(_d.loc[_d.y == 1, "p"].to_numpy()),
                        "neg": _logit(_d.loc[_d.y == 0, "p"].to_numpy())}

# Control pairs (panel A) are a separate validated set, not part of IxIDN-Matched; their p(i,j) is read from the
# full UPNA-DDAL probability matrix (stored in fig5_control_pairs.tsv). Stratum = the pair's own rarity labels.
_CP = pd.read_csv(os.path.join(FIG_DATA, "fig5_control_pairs.tsv"), sep="\t", float_precision="round_trip")
_CP = {(r.disease_a, r.disease_b): r.p for r in _CP.itertuples()}
CONTROL_PAIR_SCORES = []
for _nameA, _idA, _nameB, _idB, _role in PANEL_A_PAIRS:
    _p_ctrl = _CP[(_idA, _idB)]
    _rareA = bool(dm.loc[_idA, "is_rare"]); _rareB = bool(dm.loc[_idB, "is_rare"])
    _stratum = "rare-rare" if (_rareA and _rareB) else ("common-common" if not (_rareA or _rareB) else "rare-common")
    CONTROL_PAIR_SCORES.append(dict(name=f"{_nameA} / {_nameB}", p=_p_ctrl, logit=_logit([_p_ctrl])[0],
                                     role=_role, stratum=_stratum))
BL = [("HPO Resnik-BMA", "HPO similarity", "#7fa6cf", "s"), ("GO Resnik gene-BMA", "GO similarity", "#4c9a6a", "^")]
ST = [("rare-rare", "Rare-\nrare"), ("rare-common", "Rare-\ncommon"), ("common-common", "Common-\ncommon"), ("all", "All\npairs")]


# Panel c values as reported in the manuscript (asserted). Common-common has only 26 negatives -> drawn hollow
# (white face, edge in the series colour); its value, CI and counts are unchanged.
EXPECTED_C = {  # stratum: (auroc, ci_lo, ci_hi, HPO, GO, n_pos, n_neg)
    "rare-rare": (0.7224, 0.6716, 0.7723, 0.5753, 0.6007, 1284, 1168),
    "rare-common": (0.7580, 0.6676, 0.8346, 0.5337, 0.6090, 253, 988),
    "common-common": (0.8041, 0.6546, 0.9329, 0.5543, 0.5770, 645, 26),
    "all": (0.7645, 0.7292, 0.7993, 0.5441, 0.6446, 2182, 2182)}
HOLLOW_C = {"common-common"}   # few negatives: low-precision estimate


def draw_panel_c_allpairs(ax, panel_letter):
    R = T[T.scope == "All 4364 pairs"].set_index("stratum")
    assert [st for st, _ in ST] == list(EXPECTED_C)
    for i, (st, lab) in enumerate(ST):
        r = R.loc[st]
        got = (r.auroc, r.ci_lo, r.ci_hi, r["HPO Resnik-BMA | auroc"], r["GO Resnik gene-BMA | auroc"], r.n_pos, r.n_neg)
        assert all(abs(g - e) < 6e-5 if k < 5 else g == e for k, (g, e) in enumerate(zip(got, EXPECTED_C[st]))), (st, got)
        hol = st in HOLLOW_C
        ax.errorbar(i, r.auroc, yerr=[[r.auroc - r.ci_lo], [r.ci_hi - r.auroc]], fmt="o", color="#d9822b",
                    ms=6, capsize=3, lw=1.4, zorder=3, label="DisAtlas" if i == 0 else None,
                    mfc="white" if hol else "#d9822b", mec="#d9822b", mew=1.4)
        for k, (col_key, lab2, col, m) in enumerate(BL):
            ax.scatter(i + 0.14 + 0.07 * k, r[f"{col_key} | auroc"], marker=m, s=20, zorder=2,
                       facecolors="white" if hol else col, edgecolors=col, linewidths=1.1 if hol else 0,
                       label=lab2 if i == 0 else None)
        ax.text(i + 0.12, 0.455, f"{int(r.n_pos)}/{int(r.n_neg)}", ha="center", fontsize=9.5, color="0.35")
    ax.axhline(0.5, color="0.8", ls="--", lw=0.9)
    ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l in ST])
    ax.set_ylim(0.43, 0.98); ax.set_xlim(-0.45, 3.75)
    ax.set_ylabel("AUROC (DisAtlas: 95% CI)")
    ax.set_title("All IxIDN-Matched pairs", loc="left", fontsize=12)
    panel_label(ax, panel_letter)
    h, l = ax.get_legend_handles_labels()
    n_few = int(R.loc["common-common"].n_neg)
    h.append(Line2D([], [], color="0.35", marker="o", mfc="white", mec="0.35", mew=1.4, ms=6, ls=""))
    l.append(f"hollow = few negatives (n = {n_few})")
    ax.legend(h, l, loc="upper left", frameon=False, fontsize=10.5, bbox_to_anchor=(0.0, 1.02))


def draw_panel_d_scoredist(ax, panel_letter):
    """UPNA-DDAL score p(i,j) distributions (logit axis), true-pair vs non-pair, by rarity stratum.
    Overlaid: panel A's 5 validated control pairs as colored ticks on their own rarity stratum."""
    ST_D = [("rare-rare", "Rare-\nrare"), ("rare-common", "Rare-\ncommon"),
            ("common-common", "Common-\ncommon"), ("all", "All\npairs")]
    assert [x for x, _ in ST_D] == [x for x, _ in ST] == list(EXPECTED_C)
    OFF, WID = 0.18, 0.30
    parts_pos, parts_neg = [], []
    for i, (st, lab) in enumerate(ST_D):
        for sign, off, color in [("pos", -OFF, GOOD), ("neg", OFF, BAD)]:
            vals = SCORE_DIST[st][sign]
            vp = ax.violinplot([vals], positions=[i + off], widths=WID, showmeans=False,
                                showmedians=True, showextrema=False)
            for b in vp["bodies"]:
                b.set_facecolor(color); b.set_edgecolor(color); b.set_alpha(0.42); b.set_zorder(2)
            vp["cmedians"].set_color(color); vp["cmedians"].set_linewidth(1.6); vp["cmedians"].set_zorder(3)
            (parts_pos if sign == "pos" else parts_neg).append(vp["bodies"][0])
        npos, nneg = len(SCORE_DIST[st]["pos"]), len(SCORE_DIST[st]["neg"])
        assert (npos, nneg) == EXPECTED_C[st][5:], (st, npos, nneg)   # same strata order and counts as panel c
        ax.text(i, -11.1, f"{npos}+ / {nneg}\u2212", ha="center", fontsize=9.5, color="0.35")
    # control pairs: tick at own p(i,j), on the x-offset matching its pos/neg role, within its own
    # rarity stratum (from the stored disease labels, not inferred from name/role)
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
    ax.set_ylabel(r"UPNA-DDAL score  $p_{i,j}$  (logit scale)")
    ax.set_title("Score distribution by pair rarity", loc="left", fontsize=12)
    panel_label(ax, panel_letter)
    parts_pos[0].set_label("True pair ($y=1$)"); parts_neg[0].set_label("Non-pair ($y=0$)")
    leg_handles = [parts_pos[0], parts_neg[0],
                   Line2D([], [], color=GOOD, lw=2.6, solid_capstyle="butt", label="positive control pair"),
                   Line2D([], [], color=BAD, lw=2.6, solid_capstyle="butt", label="negative control pair")]
    ax.legend(handles=leg_handles, loc="upper center", ncol=2, fontsize=9.5,
              handlelength=1.6, columnspacing=1.1, labelspacing=0.35, handletextpad=0.5,
              framealpha=0.9, edgecolor="none", borderpad=0.6)

fig = plt.figure(figsize=(15.0, 9.6))
axA = fig.add_axes([0.02, 0.48803, 0.4533, 0.47587])
axB = fig.add_axes([0.5267, 0.48803, 0.4533, 0.47587])
axC = fig.add_axes([0.02, 0.1052, 0.43165, 0.3497])
axD = fig.add_axes([0.54165, 0.1052, 0.43165, 0.3497])

draw_panel_A_pairs(axA, PANEL_A_PAIRS, PANEL_A_SINGLES, "a", "five validated control pairs")
draw_panel_rarecommon(axB, PANEL_B_TARGETS, "b", "named clinical examples")
draw_panel_c_allpairs(axC, "c")
draw_panel_d_scoredist(axD, "d")

# write the figure straight to the manuscript FIGURES folder (no stray file in the working directory)
save(fig, "fig5_disease_map", bbox_inches="tight")
print("saved", os.path.join(HERE, "fig5_disease_map.png"))
