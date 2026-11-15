"""SI Methods figure: how the DisAtlas training graph is built from its two evidence channels and turned into a cosine score.
Every count is read from data/derived/si_graph_channels/SI_graph_channels_counts.tsv (computed from the graph build files;
the build files contain licensed edges and are not shipped). Panels a-c and the first part of panel d describe the
channel build; the rest of panels d and e give the gene-anchored training graph (rows named ga_*) that the released model
is trained on. The edge-weight ranges are read from edge_weight_range.json in the same folder. Identity checks between the
counts are re-asserted here. Panel f (provenance of the shared evidence of the four Figure 4 pairs) is drawn from the
aggregate share tables in provenance_f/ and placed under panels a-e (Lanczos resize to the width of panels a-e).
Usage: python figures/SI/si_graph_channels_schematic.py"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED  # noqa: E402
import sys, json, numpy as np, pandas as pd, matplotlib as mpl
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle, RegularPolygon, Rectangle, FancyArrowPatch
from disatlas import figstyle_si as fs
OUT = HERE; fs.apply_figure_style()
T = DERIVED / "si_graph_channels"
CT = pd.read_csv(T / "SI_graph_channels_counts.tsv", sep="\t")
C = {r.quantity: (int(r.value), "", r.how) for r in CT.itertuples()}
fmt = lambda n: f"{C[n][0]:,}"
_w = json.load(open(T / "edge_weight_range.json")); wmin, wmax = _w["w_min"], _w["w_max"]
assert C["union_edges"][0] == C["broad_edges"][0] + C["spec_edges"][0] - C["both_edges"][0]
assert C["st7_train_edges"][0] == C["train_edges"][0] == C["union_edges"][0] - C["st6_total"][0]
assert C["union_edges"][0] == C["st3_comb_edges"][0] - C["st4_ctrl_total"][0] and C["st2_proj_kept"][0] == C["spec_proj_edges"][0]
assert C["broad_edges"][0] == C["st1_broad_edges"][0] - C["st4_ctrl_broad"][0] and C["spec_edges"][0] == C["spec_proj_edges"][0] - C["st4_ctrl_spec"][0]
assert C["broad_nodes"][0] == C["st5_dis_broad"][0] and C["union_nodes"][0] == C["st7_train_dis"][0]
assert C["spec_links"][0] == sum(C[f"links_{g_}"][0] for g_ in ["gene", "lit", "path", "mech"])
assert C["ga_union_edges"][0] == C["union_edges"][0] + C["ga_added"][0]
assert C["ga_union_edges"][0] == C["ga_broad_edges"][0] + C["ga_spec_edges"][0] - C["ga_both_edges"][0]
assert C["ga_train_edges"][0] == C["ga_union_edges"][0] - C["ga_st6_total"][0] and C["ga_st6_total"][0] == C["ga_st6_pos"][0] + C["ga_st6_neg"][0]
assert C["ga_broad_edges"][0] == C["broad_edges"][0]
gwmin, gwmax = _w["ga_w_min"], _w["ga_w_max"]

# ------------------------------------------------------------------ drawing helpers
BRD, SPC, DIS, GRY, INK = "#c9971c", "#e0661b", "#2166ac", "0.55", "0.15"
MC = {"gene": "#7b3294", "lit": "#008837", "path": "#5e3c99", "mech": "#a6501b"}
MC["path"] = "#4d4d4d"
fig = mpl.pyplot.figure(figsize=(7.2, 5.75))
def panel(rect, letter, title, y0=0):
    ax = fig.add_axes(rect); ax.set_xlim(0, 1); ax.set_ylim(y0, 1); ax.axis("off")
    ax.text(0, 1.0, letter, fontsize=9, fontweight="bold", va="top", transform=ax.transAxes)
    ax.text(.028 / rect[2], 1.0, title, fontsize=7.5, va="top", transform=ax.transAxes); return ax
def asp(ax):
    b = ax.get_position(); W, H = fig.get_size_inches(); y0, y1 = ax.get_ylim(); return (b.height * H / (y1 - y0)) / (b.width * W)
def dnode(ax, x, y, t, r=None):
    r = r or (.075 if ax.get_ylim()[0] < -.4 else .06)
    ax.add_patch(Ellipse((x, y), 2 * r * asp(ax), 2 * r, fc="white", ec=DIS, lw=1.2, zorder=4)); ax.text(x, y, t, fontsize=5.5, ha="center", va="center", zorder=5)
def sq(ax, x, y, c, s=.035): ax.add_patch(Rectangle((x - s * asp(ax) / 2, y - s / 2), s * asp(ax), s, fc="white", ec=c, lw=1, zorder=4))
def hexn(ax, x, y, c, r=.033): ax.add_patch(Ellipse((x, y), 2 * r * asp(ax), 2 * r, fc="#f2f2f2", ec=c, lw=1.1, zorder=4))
def arrow(ax, a, b, **kw): ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=8, color=INK, lw=.9, **kw))
def note(ax, x, y, t, **kw): ax.text(x, y, t, **{"fontsize": 6, "va": "top", "ha": "left", "color": INK, **kw})

# a: broad channel
a = panel([.02, .575, .46, .402], "a", "Broad channel: shared genes and phenotypes", y0=-.5)
D = [(.08, .78), (.08, .58), (.08, .38)]
G = [(.33, .86), (.33, .72), (.33, .58)]; H = [(.33, .44), (.33, .30)]
for i, (x, y) in enumerate(D): dnode(a, x, y, f"D{i + 1}")
for x, y in G: sq(a, x, y, BRD)
for x, y in H: a.add_patch(Ellipse((x, y), .035 * asp(a), .035, fc=BRD, ec=BRD, lw=1, zorder=4))
for d, t in [(0, 0), (0, 1), (1, 1), (1, 2), (2, 2)]: a.plot([D[d][0], G[t][0]], [D[d][1], G[t][1]], color=BRD, lw=.8, zorder=1)
for d, t in [(0, 0), (1, 0), (1, 1), (2, 1)]: a.plot([D[d][0], H[t][0]], [D[d][1], H[t][1]], color=BRD, lw=.8, ls=(0, (3, 1.5)), zorder=1)
note(a, .37, .9, "genes (GO terms, cell-type\nexpression, protein domains)", fontsize=5.5)
note(a, .37, .46, "phenotypes\n(HPO terms)", fontsize=5.5)
arrow(a, (.6, .6), (.7, .6))
for i, (x, y) in enumerate([(.8, .78), (.92, .55), (.78, .38)]): dnode(a, x, y, f"D{i + 1}")
for p, q in [((.8, .78), (.92, .55)), ((.92, .55), (.78, .38))]: a.plot(*zip(p, q), color=BRD, lw=1.6, zorder=1)
note(a, .66, .26, "edge weight = 1", ha="left", fontsize=5.5)
for xx, ls_, t in [(.0, "-", "disease–gene link"), (.33, (0, (3, 1.5)), "disease–phenotype link")]:
    a.plot([xx, xx + .06], [.17, .17], color=BRD, lw=.8, ls=ls_); note(a, xx + .075, .17, t, fontsize=5.5, va="center")
note(a, .0, .12, f"{fmt('broad_gene_rows')} disease–gene links ({fmt('broad_gene_dis')} diseases) and\n"
     f"{fmt('broad_hpo_rows')} disease–phenotype links ({fmt('broad_hpo_dis')} diseases).\n"
     "Four similarity views (GO, phenotype, cell-type expression,\n"
     "protein domain) are rank-aggregated; the best pairs in each\n"
     "evidence class, and at least 3 per disease, are kept:\n"
     f"{fmt('st1_broad_edges')} edges on {fmt('st1_broad_dis')} diseases as built; {fmt('st4_ctrl_broad')} of\n"
     "these are later removed as FDA control pairs (see d).")

# b: specific channel
b = panel([.52, .575, .46, .402], "b", "Specific channel: shared mechanism modules", y0=-.5)
D = [(.08, .84), (.08, .66), (.08, .48), (.08, .30)]
for i, (x, y) in enumerate(D): dnode(b, x, y, f"D{i + 1}")
M = {"gene": (.42, .86), "lit": (.42, .66), "path": (.42, .48), "mech": (.42, .30)}
lab = {"gene": f"pathway modules from\ndisease genes ({fmt('mods_gene')})",
       "lit": f"literature: processes, pathways,\nmechanism entities ({fmt('mods_lit')})",
       "path": f"pathology text: anatomy, cell\ntype, phenotype ({fmt('mods_path')})",
       "mech": f"mechanism nodes ({fmt('mods_mech')}, of which\n{fmt('mods_mech_curated')} curated complement)"}
for g_, (x, y) in M.items(): hexn(b, x, y, MC[g_]); note(b, x + .05, y + .03, lab[g_], fontsize=5.5, color=MC[g_])
for d, g_, w in [(0, "gene", 1.6), (1, "gene", .8), (1, "lit", 1.2), (2, "lit", .6), (2, "path", 1.4), (3, "path", .7), (0, "mech", 1.0), (3, "mech", 1.5)]:
    b.plot([D[d][0], M[g_][0]], [D[d][1], M[g_][1]], color=MC[g_], lw=w, zorder=1)
note(b, .0, .12, f"{fmt('spec_links')} disease–module links ({fmt('spec_dis')} diseases;\n"
     f"{fmt('spec_mods')} modules, each shared by at least two diseases).\n"
     "Link weight e = number of independent source families\n"
     "+ ln(1 + dated documents). A link is kept if at least two\n"
     "source families or two dated documents support it.\n"
     "Line width = e.")

# c: projection
c = panel([.02, .222, .46, .368], "c", "Specific channel: projection to disease–disease", y0=-.3)
dnode(c, .1, .78, "D1"); dnode(c, .1, .38, "D2")
hexn(c, .4, .78, GRY); hexn(c, .4, .38, GRY)
note(c, .45, .86, "rare module\n(2 diseases)", fontsize=5.5); note(c, .45, .46, "common module\n(40 diseases)", fontsize=5.5)
for p, w in [((.4, .78), 2.4), ((.4, .38), .7)]:
    c.plot([.1, p[0]], [.78, p[1]], color=SPC, lw=w, zorder=1); c.plot([.1, p[0]], [.38, p[1]], color=SPC, lw=w, zorder=1)
arrow(c, (.68, .58), (.78, .58)); dnode(c, .88, .78, "D1"); dnode(c, .88, .38, "D2"); c.plot([.88, .88], [.735, .425], color=SPC, lw=2, zorder=1)
note(c, .91, .6, "s", fontsize=6.5, style="italic")
note(c, .0, .2, r"$x_{dm}=\sqrt{e_{dm}\,/\,(\ln(1+n_m)\,L_d)}$,  $s_{ij}=\sum_m x_{im}\,x_{jm}$" "\n"
     "$n_m$ = diseases on module m; $L_d$ = modules on disease d.\n"
     "Rare shared modules count more. Each disease keeps its 50\n"
     f"strongest partners: {fmt('st2_proj_raw')} edges; {fmt('st2_proj_kept')} on {fmt('st2_proj_dis')}\n"
     f"diseases after removing {fmt('st2_proj_drop')} that reach a disease outside the broad graph.")

# d: union and walk weights
d = panel([.52, .222, .46, .368], "d", "Combined graph and walk weights", y0=-.3)
P = {1: (.1, .8), 2: (.32, .86), 3: (.48, .62), 4: (.3, .42), 5: (.1, .5)}
for i, (x, y) in P.items(): dnode(d, x, y, f"D{i}")
for p, q, kind in [(1, 2, "b"), (2, 3, "b"), (1, 5, "b"), (5, 4, "b"), (2, 4, "s"), (3, 4, "bs"), (1, 4, "s")]:
    xy = list(zip(P[p], P[q]))
    if kind == "bs": d.plot(*xy, color=BRD, lw=3.2, zorder=1); d.plot(*xy, color=SPC, lw=1.4, zorder=2)
    else: d.plot(*xy, color=BRD if kind == "b" else SPC, lw=1.0 if kind == "b" else 2.2, zorder=1)
for y, col, t in [(.86, BRD, "broad only"), (.76, SPC, "specific only"), (.66, None, "in both")]:
    if col: d.plot([.62, .7], [y, y], color=col, lw=1.8)
    else: d.plot([.62, .7], [y, y], color=BRD, lw=3.2); d.plot([.62, .7], [y, y], color=SPC, lw=1.4)
    note(d, .72, y + .025, t, fontsize=5.5)
note(d, .0, .2, r"$w = 1[\mathrm{broad}] + 2\,f\,s\,/\,\overline{s}\;1[\mathrm{specific}]$" "\n"
     f"{fmt('st3_comb_edges')} edges; removing {fmt('st4_ctrl_total')} FDA repurposing control pairs ({fmt('st4_ctrl_broad')} broad,\n"
     f"{fmt('st4_ctrl_spec')} specific, {fmt('st4_ctrl_both')} in both) leaves {fmt('union_edges')}. Pairs sharing a\n"
     f"high-confidence causal gene add {fmt('ga_added')} specific edges and give specific\n"
     f"weight to {fmt('ga_spec_gained')} more: {fmt('ga_union_edges')} on {fmt('ga_union_nodes')} diseases ({fmt('ga_broad_edges')} broad,\n"
     f"{fmt('ga_spec_edges')} specific, {fmt('ga_both_edges')} both); w from {gwmin:.2f} to {gwmax:.1f}.\n"
     f"Balancing f: of {fmt('mech_classes')} mechanism classes, those with less walk\n"
     "weight than complement have their edges raised\n"
     f"(at most ×10; {fmt('ga_raised_edges')} edges raised); otherwise f = 1.")

# e: embedding and link score
e = panel([.02, .0, .96, .215], "e", "From graph to embedding and cosine score")
boxes = [(.0, "training graph", f"{fmt('ga_union_edges')} − {fmt('ga_st6_total')} benchmark\npairs ({fmt('ga_st6_pos')} positive,\n{fmt('ga_st6_neg')} negative) =\n{fmt('ga_train_edges')} edges"),
         (.2, "weighted random walks", "second-order, p = 1, q = 0.5;\nlength 20, window 2; w sets\nstep probabilities"),
         (.4, "disease embedding", "32 dimensions per disease;\n5 independently trained seeds"),
         (.6, "pair score", "cosine of the two disease\nvectors in each seed;\nno trained parameters"),
         (.8, "cosine score", "mean of the 5 per-seed\ncosines; used to rank\ncandidate disease pairs")]
for x, h, t in boxes:
    e.add_patch(FancyBboxPatch((x + .01, .16), .155, .56, boxstyle="round,pad=0.01,rounding_size=0.02", fc="#f4f6fa", ec=DIS, lw=.9))
    e.text(x + .0875, .68, h, fontsize=6.1, fontweight="bold", ha="center", va="top")
    e.text(x + .0875, .5, t, fontsize=5.4, ha="center", va="top", color=INK)
for x, _, _ in boxes[:-1]: arrow(e, (x + .178, .44), (x + .202, .44))
fig.text(.5, .0, "", fontsize=1)
N = {n: v for n, (v, _, _) in C.items()}; N.update(w_min=float(wmin), w_max=float(wmax), ga_w_min=float(gwmin), ga_w_max=float(gwmax))
fs.save(fig, OUT, "SI_graph_channels_schematic_ae"); N["overlaps"] = fs.overlap_check(fig)
mpl.pyplot.close(fig)

# f: provenance of the shared evidence of the four Figure 4 pairs (aggregate shares; drawn as a separate 180 mm panel)
PF = T / "provenance_f"
fs.apply_figure_style(sizes=(7, 6, 6))
P4 = ["aHUS-PNH", "NMOSD-PNH", "Crohn-UC", "RA-polymyositis"]
MP = [("D–G–Pathway–G–D", "Genes in a shared pathway"), ("D–G–GeneSet–G–D", "Genes in a shared named gene set"),
      ("D–G–D shared rare gene", "Shared rare-disease gene"),
      ("D–Mechanism node–D", "Shared mechanism node"), ("D–Mechanism entity–D", "Shared literature mechanism entity"),
      ("D–Pathway–D", "Shared literature pathway"), ("D–GO process–D", "Shared GO process"),
      ("D–Phenotype–D", "Shared phenotype"), ("D–Cell type–D", "Shared cell type"), ("D–Anatomy–D", "Shared anatomy"),
      ("D–D broad similarity", "Broad similarity")]
FM = ["Curated databases", "Expression atlas", "Literature text mining", "Pathology text mining", "Broad similarity"]
A = pd.read_csv(PF / "metapath_share_matrix.tsv", sep="\t", index_col=0)
B = pd.read_csv(PF / "source_family_share_matrix.tsv", sep="\t", index_col=0)
FT = pd.read_csv(PF / "pair_source_family.tsv", sep="\t")
assert list(A.index) == [m for m, _ in MP] and list(B.index) == FM and list(A.columns) == P4 == list(B.columns)
assert np.allclose(A.sum(), 1) and np.allclose(B.sum(), 1)
assert np.allclose(FT.pivot(index="family", columns="pair", values="share").reindex(index=FM, columns=P4).fillna(0).values, B.values)
HB = FT.groupby("pair").hub_item_share.sum().reindex(P4)
mm = 25.4; cmap = mpl.colormaps["Blues"]; norm = mpl.colors.Normalize(0, 0.40)
PL2 = ["aHUS–\nPNH", "NMOSD–\nPNH", "Crohn’s–\nUC", "RA–poly-\nmyositis"]
def hm2(ax, M, ylab):
    ax.imshow(M.values, cmap=cmap, norm=norm, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]; t = "–" if v == 0 else ("<1" if v < 0.005 else f"{100*v:.0f}")
            ax.text(j, i, t, ha="center", va="center", fontsize=6, color="white" if norm(v) > 0.55 else "#222")
    ax.set_xticks(range(4)); ax.set_xticklabels(PL2); ax.xaxis.tick_top()
    ax.set_yticks(range(len(ylab))); ax.set_yticklabels(ylab); ax.tick_params(length=0, pad=2)
    for s_ in ax.spines.values(): s_.set_visible(False)
ff = mpl.pyplot.figure(figsize=(180/mm, 64/mm), dpi=300)
a1 = ff.add_axes([0.215, 0.03, 0.235, 0.78]); a2 = ff.add_axes([0.62, 0.03, 0.235, 0.78]); cax = ff.add_axes([0.875, 0.12, 0.012, 0.6])
hm2(a1, A, [l for _, l in MP]); a1.axhline(len(MP) - 1.5, color="white", lw=2)
hm2(a2, B, FM); a2.axhline(3.5, color="white", lw=2)
for j, p_ in enumerate(P4): a2.text(j, 5.0, f"{100*HB[p_]:.0f}", ha="center", va="center", fontsize=6, color="#555")
a2.set_yticks(list(range(5)) + [5.0]); a2.set_yticklabels(FM + ["Of which hub modules"]); a2.get_yticklabels()[-1].set_color("#555"); a2.set_ylim(5.4, -0.5)
cb = ff.colorbar(mpl.cm.ScalarMappable(norm, cmap), cax=cax, extend="max"); cb.set_ticks([0, .1, .2, .3, .4]); cb.set_ticklabels(["0", "10", "20", "30", "40"])
cb.set_label("Share of shared walk mass (%)", fontsize=6); cb.outline.set_visible(False)
a1.set_title("By shared evidence type", loc="left", fontsize=7, pad=20); a2.set_title("By source family", loc="left", fontsize=7, pad=20)
ff.canvas.draw(); rr = ff.canvas.get_renderer(); bb = a1.title.get_window_extent(rr).transformed(ff.transFigure.inverted())
ff.text(0.008, bb.y1, "f", fontsize=9, fontweight="bold", va="top")
ff.savefig(f"{OUT}/SI_graph_channels_schematic_f.png", dpi=300); ff.savefig(f"{OUT}/SI_graph_channels_schematic_f.pdf", metadata={"CreationDate": None})
N["overlaps_f"] = fs.overlap_check(ff); mpl.pyplot.close(ff)

# assemble: panels a-e on top, 42 px white gap, panel f resized to the same width
from PIL import Image
TOP = Image.open(f"{OUT}/SI_graph_channels_schematic_ae.png").convert("RGB")
PFI = Image.open(f"{OUT}/SI_graph_channels_schematic_f.png").convert("RGB")
GAP = 42; fh = round(PFI.size[1] * TOP.size[0] / PFI.size[0])
CMP = Image.new("RGB", (TOP.size[0], TOP.size[1] + GAP + fh), "white")
CMP.paste(TOP, (0, 0)); CMP.paste(PFI.resize((TOP.size[0], fh), Image.LANCZOS), (0, TOP.size[1] + GAP))
CMP.save(f"{OUT}/SI_graph_channels_schematic.png", dpi=(300, 300))
_pf = mpl.pyplot.figure(figsize=(CMP.size[0] / 300, CMP.size[1] / 300), dpi=300); _pf.figimage(np.asarray(CMP), 0, 0)
_pf.savefig(f"{OUT}/SI_graph_channels_schematic.pdf", dpi=300, metadata={"CreationDate": None}); mpl.pyplot.close(_pf)  # raster PDF of the assembled figure; vector panels in _ae.pdf and _f.pdf
import os; N["png_bytes"] = os.path.getsize(f"{OUT}/SI_graph_channels_schematic.png")
json.dump(N, open(f"{OUT}/SI_graph_channels_schematic_numbers.json", "w"), indent=1, default=str); print(json.dumps(N, default=str))
