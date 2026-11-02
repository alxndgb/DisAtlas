"""Main Figure 3: the IxIDN-Matched benchmark (four panels).
  a  how IxIDN is built: a real four-disease, four-drug excerpt of the disease-drug trial graph (black lines)
     projected onto diseases (gold IxIDN links); the three links that are IxIDN+ pairs of the benchmark in green;
  b  the construction steps from IxIDN to IxIDN-Matched, with the pair count at every stage (the full funnel is
     Supplementary Figure SI_ixidn_matched_funnel);
  c  how an IxIDN- pair is chosen, for the real matched pair E01438 (IxIDN+) / E01435 (IxIDN-, same matching cell);
  d  single-property AUROCs of the matched covariates, with chance (0.5) and DisAtlas (0.832, the full model on IxIDN-Matched) for reference.

Values are drawn as labels; check_sources() re-reads the repository tables they come from and stops if any differs.
Inputs (data/derived/fig3_ixidn/ unless marked [Z] = data/zenodo/figures_main/):
  fig3a_evidence.csv, dd_trial_graph_excerpt.tsv              (panel a disease-drug lines and trial counts)
  [Z] ixidn_graph_edges.parquet                               (panel a gold links; panel c shortest paths)
  funnel_counts.tsv + stage1_ixidn_pos.tsv.gz, stage2_universe_diseases.tsv.gz, stage3_ixidn_pos.tsv.gz,
  [Z] stage3_ixidn_neg.tsv.gz                                 (panel b: 1,379 / 43,112, 985, 27,415)
  study_intensity_profile.tsv                                 (panel b: 155,375)
  [Z] stage3_candidates_degree.parquet                        (panel b: degree-sum AUROC 0.947 on the stage-3 candidates)
  matched_candidate_verdicts.tsv, ixidn_matched_pairs.tsv     (panel a green links; panel b 4,302 -> 2,182 + 2,182, 938; panel c)
  fig3c_distances.csv                                         (panel c shortest paths)
  shortcut_diagnostics.csv                                    (panel d single-property AUROCs)
Writes ixidn_benchmark.png/.pdf and ixidn_benchmark_values.tsv next to this script.

    python figures/main/fig3_ixidn_benchmark.py [out.png]
"""
import os, sys, csv, json, collections
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
X = os.path.join(ROOT, "data", "derived", "fig3_ixidn")
Z = os.path.join(ROOT, "data", "zenodo", "figures_main")
from disatlas.figstyle import apply_style, panel_label, LEGEND_SIZE, ANNOT_SIZE, BASE_SIZE, TITLE_SIZE, DPI
from disatlas.parity import write_values
apply_style()

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon, FancyBboxPatch, Circle
from matplotlib.lines import Line2D

# ---------------------------------------------------------------- drawn values (checked in check_sources)
DISEASES = ["Lewy body dementia", "Parkinson disease", "Alzheimer disease", "alcohol abuse"]
DRUGS = ["pimavanserin", "creatine", "trazodone", "prazosin"]
DD_LINES = [("Lewy body dementia", "pimavanserin"), ("Parkinson disease", "pimavanserin"), ("Parkinson disease", "creatine"),
            ("Parkinson disease", "trazodone"), ("Alzheimer disease", "pimavanserin"), ("Alzheimer disease", "creatine"),
            ("Alzheimer disease", "trazodone"), ("Alzheimer disease", "prazosin"), ("alcohol abuse", "trazodone"),
            ("alcohol abuse", "prazosin")]
# projected IxIDN links; third field = the IxIDN+ benchmark pair drawn in green (None = gold)
LINKS = [("Lewy body dementia", "Alzheimer disease", "E01438"), ("Lewy body dementia", "Parkinson disease", "E01422"),
         ("Parkinson disease", "alcohol abuse", "E03971"), ("Parkinson disease", "Alzheimer disease", None),
         ("Alzheimer disease", "alcohol abuse", None)]
FUNNEL = dict(ixidn_diseases=1379, ixidn_links=43112, universe_diseases=985, cand_pos=27415, cand_neg=155375,
              dsum_auroc="0.947", matched_pos=4302, bench_pos=2182, bench_neg=2182, bench_diseases=938)
PAIR_C = dict(pos="E01438", neg="E01435", cell="0|1|1|-1|1|3|1", drug="pimavanserin",
              u="Alzheimer disease", v="Lewy body dementia", u2="stroke disorder",
              w="epilepsy with myoclonic atonic seizures", dist={"u": 1, "v": 2, "w": 3})
SHORTCUTS = [("MONDO similarity", "mondo_simgic", 0.5418), ("degree sum", "dsum", 0.5282),
             ("name similarity", "name_tfidf_cos", 0.5062), ("gene count (sum, log)", "lgp", 0.5043),
             ("gene count (smaller)", "gmin", 0.5011)]
import pandas as _pd
from sklearn.metrics import roc_auc_score as _auc
_R = _pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "derived", "fig6_comparison", "per_run_scores.tsv"), sep="\t")
DISATLAS_AUROC = round(float(_auc(_R.y, _R["DisAtlas"])), 3)  # IxIDN-Matched full-model AUROC from the shipped scores

# ---------------------------------------------------------------- colours and sizes
PUR = ("#e9e0f3", "#6a4a9c")      # disease nodes
DRG = ("#e1f1f8", "#3d87a8")      # drug nodes
GOLD = "#c9930e"                  # IxIDN link
POS = "#3f7f52"                   # IxIDN+ (panels a, b, c)
NEG = "#b03a3a"                   # IxIDN- (panels b, c)
NEG_DARK = "#7a2b2b"
BLK = "#222222"
BOXF, BOXE = "#eef0f5", "#5a6080"
SUB = "#444444"
FS, FL, FA = BASE_SIZE, LEGEND_SIZE, ANNOT_SIZE    # body / node and legend text / annotations

fig = plt.figure(figsize=(15.3, 13.4))
ax_a = fig.add_axes([0.035, 0.640, 0.505, 0.300])
ax_c = fig.add_axes([0.035, 0.035, 0.505, 0.545])
ax_b = fig.add_axes([0.555, 0.335, 0.430, 0.605])
ax_d = fig.add_axes([0.700, 0.055, 0.285, 0.200])
LABEL_X = 0.012                                    # figure x of the panel letters in the left column
LABEL_X_R = 0.535                                  # ... and in the right column


def head(ax, letter, title, fx, fy):
    """House panel letter at figure position (fx, fy) and a short plain title beside it."""
    p = ax.get_position()
    panel_label(ax, letter, xy=((fx - p.x0) / p.width, (fy - p.y0) / p.height))
    fig.text(fx + 0.020, fy, title, fontsize=TITLE_SIZE, va="bottom", ha="left")


# ---------------------------------------------------------------- panel a
def draw_a(ax):
    ax.set_xlim(0, 100); ax.set_ylim(22, 95); ax.axis("off")
    DX, GX = 27, 40
    dy = dict(zip(DISEASES, (84, 69, 54, 39))); gy = dict(zip(DRUGS, (84, 69, 54, 39)))
    for d, g in DD_LINES:
        ax.plot([DX, GX], [dy[d], gy[g]], color=BLK, lw=1.0, zorder=1)
    for d, y in dy.items():
        ax.scatter(DX, y, s=110, c=PUR[0], edgecolors=PUR[1], zorder=3)
        ax.text(DX - 2.2, y, d, ha="right", va="center", fontsize=FL)
    for g, y in gy.items():
        ax.scatter(GX, y, s=110, c=DRG[0], edgecolors=DRG[1], zorder=3)
        ax.text(GX + 2.2, y, g, ha="left", va="center", fontsize=FL)
    ax.text(DX, 90.5, "Disease", ha="center", va="bottom", fontweight="bold", fontsize=FL)
    ax.text(GX, 90.5, "Drug", ha="center", va="bottom", fontweight="bold", fontsize=FL)
    ax.text(22, 31, "black line = drug trialled for the disease", ha="center", va="top", fontsize=FA,
            style="italic", color=SUB)
    ax.add_patch(FancyArrowPatch((50.5, 61.5), (56.5, 61.5), arrowstyle="simple,head_width=7,head_length=6", color="#888888"))
    P = {"Lewy body dementia": (80, 84), "Parkinson disease": (71, 61.5), "Alzheimer disease": (89, 61.5),
         "alcohol abuse": (80, 42)}
    LAB = {"Lewy body dementia": ("Lewy body dementia", 80, 88.0, "center", "bottom"),
           "Parkinson disease": ("Parkinson", 68.5, 61.5, "right", "center"),
           "Alzheimer disease": ("Alzheimer", 91.5, 61.5, "left", "center"),
           "alcohol abuse": ("alcohol abuse", 83, 42, "left", "center")}
    for a, b, uid in LINKS:
        (x1, y1), (x2, y2) = P[a], P[b]
        ax.plot([x1, x2], [y1, y2], color=POS if uid else GOLD, lw=3.0 if uid else 2.2, zorder=1)
    for k, (x, y) in P.items():
        ax.scatter(x, y, s=110, c=PUR[0], edgecolors=PUR[1], zorder=3)
        t, tx, ty, ha, va = LAB[k]
        ax.text(tx, ty, t, ha=ha, va=va, fontsize=FL)
    ly = 32
    ax.plot([64, 69], [ly, ly], color=GOLD, lw=2.2); ax.text(70, ly, "IxIDN link (shared trial drug)", va="center", fontsize=FA)
    ax.plot([64, 69], [ly - 5.5, ly - 5.5], color=POS, lw=3.0)
    ax.text(70, ly - 5.5, "IxIDN\u207a pair (benchmark positive)", va="center", fontsize=FA)


# ---------------------------------------------------------------- panel b
def trap(ax, xt0, xt1, xb0, xb1, y0, y1):
    ax.add_patch(Polygon([(xt0, y1), (xt1, y1), (xb1, y0), (xb0, y0)], closed=True, fc=BOXF, ec=BOXE, lw=1.0))


def draw_b(ax):
    f = FUNNEL
    ax.set_xlim(0, 100); ax.set_ylim(-7, 100); ax.axis("off")
    trap(ax, 2, 98, 7, 93, 88, 100)
    ax.text(50, 96.2, "IxIDN", ha="center", va="center", fontsize=FS, fontweight="bold")
    ax.text(50, 91.4, f"{f['ixidn_diseases']:,} diseases, {f['ixidn_links']:,} IxIDN\u207a links", ha="center", va="center", fontsize=FL)
    trap(ax, 7, 93, 12, 88, 75, 87)
    ax.text(50, 83.2, "Evaluation universe", ha="center", va="center", fontsize=FS, fontweight="bold")
    ax.text(50, 78.4, f"{f['universe_diseases']:,} diseases", ha="center", va="center", fontsize=FL)
    ax.text(50, 68.5, f"{f['cand_pos']:,} IxIDN\u207a and {f['cand_neg']:,} IxIDN\u207b candidate pairs", ha="center", va="center",
            fontsize=FA + 0.5, style="italic", color=NEG_DARK)
    trap(ax, 12, 88, 18, 82, 29, 62)
    ax.text(50, 59.6, "Match on 7 keys", ha="center", va="top", fontsize=FS, fontweight="bold")
    keys = ["degree decile, lower-degree disease", "degree decile, higher-degree disease", "shared organ system",
            "shared Orphanet clinical area", "name similarity (quartile)", "MONDO similarity (quartile)",
            "gene-coverage category"]
    for i, k in enumerate(keys):
        ax.text(50, 54.6 - 3.2 * i, f"{i + 1}  {k}", ha="center", va="center", fontsize=FL)
    ax.text(50, 31.4, f"{f['matched_pos']:,} IxIDN\u207a pairs with a matching IxIDN\u207b", ha="center", va="center",
            fontsize=FL, fontweight="bold")
    trap(ax, 20, 80, 26, 74, 16, 26.5)
    ax.text(50, 24.4, "Language-model QC", ha="center", va="top", fontsize=FS, fontweight="bold")
    ax.text(50, 17.2, f"weak shared-drug evidence dropped\n{f['bench_pos']:,} IxIDN\u207a pairs kept", ha="center", va="bottom", fontsize=FL, linespacing=1.25)
    for x0, col, fc, lab, n in [(20, POS, "#e4f1e6", "IxIDN\u207a", f["bench_pos"]), (52, NEG, "#f6e3e3", "IxIDN\u207b", f["bench_neg"])]:
        ax.add_patch(FancyBboxPatch((x0, 1.5), 28, 11, boxstyle="round,pad=0,rounding_size=0.6", fc=fc, ec=col, lw=1.0))
        ax.text(x0 + 14, 10.2, lab, ha="center", va="center", fontsize=FS + 1, fontweight="bold", color=col)
        ax.text(x0 + 14, 4.8, f"{n:,} pairs", ha="center", va="center", fontsize=FL + 0.5)
    ax.text(50, -3.6, f"IxIDN-Matched: {f['bench_pos'] + f['bench_neg']:,} pairs over {f['bench_diseases']:,} diseases",
            ha="center", va="center", fontsize=FS, fontweight="bold")


# ---------------------------------------------------------------- panel c
def draw_c(ax):
    p = PAIR_C
    ax.set_xlim(-46, 46); ax.set_ylim(-40, 62); ax.set_aspect("equal"); ax.axis("off")
    wrap = {p["u"]: "Alzheimer\ndisease", p["v"]: "Lewy body\ndementia", p["u2"]: "stroke\ndisorder",
            p["w"]: "epilepsy with myoclonic\natonic seizures"}
    pairs = [("u", "v", p["u"], p["v"], POS, "IxIDN\u207a pair", f"co-trialled ({p['drug']})", "-"),
             ("u\u2032", "w", p["u2"], p["w"], NEG, "IxIDN\u207b pair", "\u22653 steps apart", (0, (3, 2)))]
    for x0, (a, b, na, nb, col, lab, sub, ls) in zip((-44, 8), pairs):
        y = 48
        ax.plot([x0 + 3, x0 + 33], [y, y], color=col, lw=1.8, ls=ls, zorder=1)
        ax.scatter([x0 + 3, x0 + 33], [y, y], s=520, c=[PUR[0]] * 2, edgecolors=PUR[1], zorder=3)
        ax.text(x0 + 3, y, a, ha="center", va="center", fontsize=FS, fontweight="bold", zorder=4)
        ax.text(x0 + 33, y, b, ha="center", va="center", fontsize=FS, fontweight="bold", zorder=4)
        ax.text(x0 + 18, y + 7.2, lab, ha="center", va="bottom", fontsize=FS, color=col, fontweight="bold")
        ax.text(x0 + 18, y + 2.8, sub, ha="center", va="bottom", fontsize=FA, style="italic", color=col)
        ax.text(x0 + 3, y - 3.6, wrap[na], ha="center", va="top", fontsize=FA, linespacing=1.1)
        ax.text(x0 + 33, y - 3.6, wrap[nb], ha="center", va="top", fontsize=FA, linespacing=1.1)
    R1, R2, R3 = 10, 21, 32; cy = 2
    for r, lab in [(R1, "1 step"), (R2, "2 steps"), (R3, "\u22653 steps")]:
        ax.add_patch(Circle((0, cy), r, fill=False, ls="--", ec="#bbbbbb", lw=0.9))
        ax.text(0, cy + r + 0.8, lab, ha="center", va="bottom", fontsize=FA, color="#777777")
    pt = lambda r, a: (r * np.cos(np.radians(a)), cy + r * np.sin(np.radians(a)))
    n1 = [pt(R1, a) for a in (20, 140, 285)]; n2 = [pt(R2, a) for a in (0, 60, 120, 200, 300)]
    n3 = [pt(R3, a) for a in (-25, 35, 155, 215, 255)]
    O = (0, cy)
    ed = [(O, n1[0]), (O, n1[1]), (O, n1[2]), (n1[0], n2[0]), (n1[0], n2[1]), (n1[1], n2[2]), (n1[1], n2[3]),
          (n1[2], n2[4]), (n2[0], n3[0]), (n2[3], n3[2]), (n2[4], n3[4])]
    for a, b in ed:
        ax.plot([a[0], b[0]], [a[1], b[1]], color=GOLD, lw=1.3, zorder=1)
    for q in n1: ax.scatter(*q, s=100, c="#e3f1e6", edgecolors=POS, zorder=3)
    for q in n2: ax.scatter(*q, s=100, c="#eeeeee", edgecolors="#999999", zorder=3)
    for q in (n3[0], n3[2], n3[4]):
        ax.plot([0, q[0]], [cy, q[1]], color=NEG, ls=(0, (3, 2)), lw=1.4, zorder=2)
        ax.scatter(*q, s=120, c="#f6e3e3", edgecolors=NEG, zorder=3)
    wx, wy = n3[0]
    ax.scatter(wx, wy, s=520, c="#f6e3e3", edgecolors=NEG, lw=2, zorder=5)
    ax.text(wx, wy, "w", ha="center", va="center", fontsize=FS, fontweight="bold", color=NEG, zorder=6)
    ax.scatter(*O, s=520, c=PUR[0], edgecolors=PUR[1], zorder=4)
    ax.text(0, cy, "u\u2032", ha="center", va="center", fontsize=FS, fontweight="bold", zorder=5)
    U, V = n1[0], n2[1]            # u one step and v two steps from u' (PAIR_C["dist"])
    ax.plot([U[0], V[0]], [U[1], V[1]], color=POS, lw=3.0, zorder=2)
    for q, t in ((U, "u"), (V, "v")):
        ax.scatter(*q, s=520, c=PUR[0], edgecolors=POS, lw=2, zorder=5)
        ax.text(*q, t, ha="center", va="center", fontsize=FS, fontweight="bold", zorder=6)
    hs = [Line2D([], [], color=GOLD, lw=1.4, label="IxIDN link"),
          Line2D([], [], marker="o", ls="", ms=8, mfc="#eeeeee", mec="#999999", label="2 steps: not used"),
          Line2D([], [], marker="o", ls="", ms=8, mfc="#e3f1e6", mec=POS, label="1 step from u\u2032"),
          Line2D([], [], color=NEG, ls=(0, (3, 2)), marker="o", ms=8, mfc="#f6e3e3", mec=NEG, label="IxIDN\u207b candidate")]
    ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.04), fontsize=FL, frameon=False, ncol=2,
              handlelength=2.6, handletextpad=0.6, columnspacing=2.4, labelspacing=0.7)


# ---------------------------------------------------------------- panel d
def draw_d(ax):
    feat = [s[0] for s in SHORTCUTS]; au = [s[2] for s in SHORTCUTS]
    yy = np.arange(len(feat))[::-1]
    ax.barh(yy, [a - 0.5 for a in au], left=0.5, height=0.6, color="#9a9a9a")
    for y, a in zip(yy, au):
        ax.text(a + 0.004, y, "%.3f" % a, va="center", fontsize=FL)
    ax.axvline(0.5, color="#555555", lw=0.9, ls="--")
    ax.axvline(DISATLAS_AUROC, color="#3a64a0", lw=1.4)
    ax.text(DISATLAS_AUROC - 0.006, yy[0] + 0.1, "DisAtlas\n%.3f" % DISATLAS_AUROC, color="#3a64a0", fontsize=FL,
            ha="right", va="center")
    ax.set_yticks(yy); ax.set_yticklabels(feat, fontsize=FL); ax.set_xlim(0.45, 0.85 if DISATLAS_AUROC > 0.79 else 0.80)
    ax.set_xlabel("AUROC for IxIDN\u207a vs IxIDN\u207b (0.5 = chance)", fontsize=FL)
    ax.tick_params(labelsize=FL)
    ax.text(-0.509, 1.05, "organ system, Orphanet area and gene coverage: identical by construction", transform=ax.transAxes,
            ha="left", va="bottom", fontsize=FA, color=SUB)


def build():
    draw_a(ax_a); draw_b(ax_b); draw_c(ax_c); draw_d(ax_d)
    head(ax_a, "a", "IxIDN: trial graph projected onto diseases", LABEL_X, 0.958)
    head(ax_b, "b", "From IxIDN to IxIDN-Matched", LABEL_X_R, 0.958)
    head(ax_c, "c", "How an IxIDN\u207b pair is chosen", LABEL_X, 0.585)
    head(ax_d, "d", "Matched properties alone do not separate positives from negatives", LABEL_X_R, 0.300)


# ---------------------------------------------------------------- source checks
def _auroc(y, s):
    """Mann-Whitney AUROC with average ranks for ties."""
    import pandas as pd
    y = np.asarray(y); r = pd.Series(np.asarray(s, float)).rank(method="average").values
    n1 = int(y.sum()); n0 = len(y) - n1
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def check_sources():
    """Re-read the repository tables behind panels a-d and stop if a drawn value differs. Returns the recomputed AUROC."""
    import pandas as pd, networkx as nx
    # ---- panel a: disease-drug lines against the trial summary, links against IxIDN and the benchmark
    ev = pd.read_csv(os.path.join(X, "fig3a_evidence.csv"))
    dd_ev = ev[ev.kind == "disease-drug"]
    assert {(r.a, r.b.lower()) for r in dd_ev.itertuples()} == set(DD_LINES) and len(dd_ev) == len(DD_LINES)
    dis_id = dict(zip(dd_ev.a, dd_ev.a_id)); drug_id = {b.lower(): i for b, i in zip(dd_ev.b, dd_ev.b_id)}
    assert set(dis_id) == set(DISEASES) and set(drug_id) == set(DRUGS)
    dd = pd.read_csv(os.path.join(X, "dd_trial_graph_excerpt.tsv"), sep="\t")
    sub = dd[dd.disease_id.isin(dis_id.values()) & dd.compound_id.isin(drug_id.values())]
    got = {(a, b): int(n) for a, b, n in zip(sub.disease_id, sub.compound_id, sub.num_trials)}
    for r in dd_ev.itertuples():
        k = (r.a_id, r.b_id)
        assert k in got and got[k] == int(r.n_trials), ("disease-drug line", r.a, r.b, got.get(k), r.n_trials)
    assert len(got) == len(DD_LINES), ("trial-graph lines among the four diseases and drugs not drawn",
                                       sorted(set(got) - {(r.a_id, r.b_id) for r in dd_ev.itertuples()}))
    e = pd.read_parquet(os.path.join(Z, "ixidn_graph_edges.parquet"))
    eix = {frozenset((a, b)): set(json.loads(s)) for a, b, s in zip(e.disease_1, e.disease_2, e.shared_drugs)}
    drugs_of = collections.defaultdict(set)
    for d, g in DD_LINES:
        drugs_of[d].add(drug_id[g])
    expected = {frozenset((a, b)) for i, a in enumerate(DISEASES) for b in DISEASES[i + 1:] if drugs_of[a] & drugs_of[b]}
    assert expected == {frozenset((a, b)) for a, b, _ in LINKS}, "drawn links are not the projection of the drawn lines"
    for a, b, _ in LINKS:
        k = frozenset((dis_id[a], dis_id[b]))
        assert k in eix and drugs_of[a] & drugs_of[b] <= eix[k], ("gold link not an IxIDN edge with the drugs shown", a, b)
    M = pd.read_csv(os.path.join(X, "ixidn_matched_pairs.tsv"), sep="\t", dtype={"cellA": str})
    by = M.set_index("pair_uid")
    bench_pos = set(map(frozenset, zip(M[M.label == 1].a, M[M.label == 1].b)))
    for a, b, uid in LINKS:
        k = frozenset((dis_id[a], dis_id[b]))
        if uid:
            r = by.loc[uid]
            assert int(r.label) == 1 and frozenset((r.a, r.b)) == k, ("green link", uid)
        else:
            assert k not in bench_pos, ("gold link is an IxIDN+ benchmark pair", a, b)
    # ---- panel b: funnel
    f = FUNNEL
    fc = list(csv.DictReader(open(os.path.join(X, "funnel_counts.tsv")), delimiter="\t"))
    fcv = lambda st, cl: [r for r in fc if r["stage"] == st and r["class"] == cl]
    s1 = fcv("1_graph", "IxIDN+"); assert len(s1) == 1 and int(s1[0]["pairs"]) == f["ixidn_links"], s1
    g1 = pd.read_csv(os.path.join(X, s1[0]["file"]), sep="\t")
    assert len(g1) == f["ixidn_links"] and len(set(g1.disease_a) | set(g1.disease_b)) == f["ixidn_diseases"]
    s2 = fcv("2_universe_coverage", "diseases"); assert int(s2[0]["pairs"]) == f["universe_diseases"], s2
    s3 = fcv("3_pinned_universe", "IxIDN+"); assert int(s3[0]["pairs"]) == f["cand_pos"], s3
    si = {r["stratum"]: r for r in csv.DictReader(open(os.path.join(X, "study_intensity_profile.tsv")), delimiter="\t")}
    assert int(si["universe negatives"]["n"]) == f["cand_neg"] and int(si["universe positives"]["n"]) == f["cand_pos"]
    F = pd.read_parquet(os.path.join(Z, "stage3_candidates_degree.parquet"))
    assert int((F.label == 1).sum()) == f["cand_pos"] and int((F.label == 0).sum()) == f["cand_neg"]
    pos3 = pd.read_csv(os.path.join(X, s3[0]["file"]), sep="\t", usecols=["a", "b"])
    neg3 = pd.read_csv(os.path.join(Z, fcv("3_pinned_universe", "IxIDN-")[0]["file"]), sep="\t", usecols=["a", "b"])
    key = lambda t: set(map(frozenset, zip(t.a, t.b)))
    assert key(F[F.label == 1]) == key(pos3) and key(F[F.label == 0]) <= key(neg3), "candidates are not stage-3 pairs"
    au = _auroc(F.label.values, (F.deg_a + F.deg_b).values)
    assert "%.3f" % au == f["dsum_auroc"], ("degree-sum AUROC", au)
    with open(os.path.join(X, "matched_candidate_verdicts.tsv")) as fh:
        cand = list(csv.DictReader(fh, delimiter="\t"))
    v = collections.Counter(r["verdict"] for r in cand)
    assert len(cand) == f["matched_pos"] and v["STRONG"] + v["PLAUSIBLE"] == f["bench_pos"] and v["STRONG"] == 1285 and v["PLAUSIBLE"] == 897, (len(cand), v)   # legend b
    lab = M.label.value_counts().to_dict()
    assert lab == {1: f["bench_pos"], 0: f["bench_neg"]} and len(set(M.a) | set(M.b)) == f["bench_diseases"], lab
    # ---- panel c: the matched pair and its shortest paths on the evaluation-universe IxIDN graph
    p = PAIR_C
    P, N = by.loc[p["pos"]], by.loc[p["neg"]]
    assert int(P.label) == 1 and int(N.label) == 0 and P.cellA == N.cellA == p["cell"], (P.cellA, N.cellA)
    assert (P.disease_a_name, P.disease_b_name, N.disease_a_name, N.disease_b_name) == (p["u"], p["v"], p["u2"], p["w"])
    assert drug_id[p["drug"]] in eix[frozenset((P.a, P.b))]
    uni = pd.read_csv(os.path.join(X, s2[0]["file"]), sep="\t").disease_id
    G = nx.Graph(); G.add_edges_from(zip(e.disease_1, e.disease_2)); Gu = G.subgraph(set(uni))
    dist = {"u": nx.shortest_path_length(Gu, N.a, P.a), "v": nx.shortest_path_length(Gu, N.a, P.b),
            "w": nx.shortest_path_length(Gu, N.a, N.b)}
    assert dist == p["dist"], dist
    cd = pd.read_csv(os.path.join(X, "fig3c_distances.csv"))
    row = cd[(cd["pos"] == p["pos"]) & (cd["neg"] == p["neg"])].iloc[0]
    assert (row["n1-u"], row["n1-v"], row["n1-n2"]) == (dist["u"], dist["v"], dist["w"]), row
    # ---- panel d
    diag = {r["feature"]: float(r["expanded_auroc"])
            for r in csv.DictReader(open(os.path.join(X, "shortcut_diagnostics.csv")))}
    for _, k, a in SHORTCUTS:
        assert abs(diag[k] - a) < 1e-9, (k, diag[k], a)
    print(f"check_sources: panel a {len(DD_LINES)} lines, {len(LINKS)} links ({sum(1 for l in LINKS if l[2])} IxIDN+); "
          f"panel b degree-sum AUROC {au:.6f} on {f['cand_pos']:,} + {f['cand_neg']:,} candidate pairs; "
          f"panel c distances {dist}; panel d {len(SHORTCUTS)} AUROCs: OK")
    return au


if __name__ == "__main__":
    au = check_sources()
    build()
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "ixidn_benchmark.png")
    fig.savefig(out, dpi=DPI, facecolor="white"); fig.savefig(out[:-4] + ".pdf", facecolor="white")
    f = FUNNEL; V = dict(f)
    V["dsum_auroc_recomputed"] = round(au, 6); V["disatlas_auroc_reference"] = DISATLAS_AUROC
    V.update({f"shortcut_auroc|{k}": a for _, k, a in SHORTCUTS}); V["verdict_strong"] = 1285; V["verdict_plausible"] = 897
    write_values(out[:-4] + "_values.tsv", V)
    print("wrote", out)
