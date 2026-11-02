"""Main Figure 6: cross-system comparison on the 4,364 blinded IxIDN-Matched pairs (house style).

Panels:
  a ROC curves (all 4,364 pairs; multi-run methods = mean over runs; unscored pairs ranked last)
  b ranking metrics (AUPRC, F1 at k = #positives, precision top 10%, NPV bottom 10%; 95% CI disease-cluster bootstrap;
    FIG6_CI=pair switches to a pair bootstrap)
  c run-to-run Spearman agreement of pair scores
  d Jaccard overlap of top-10% pairs
Two DisAtlas arms: full model and node-disjoint (benchmark diseases held out of training; 5 folds x 5 seeds).
Comparator scores are the stored outputs of each system (Table 1 protocols); no language-model call is made here.
Inputs: data/derived/fig6_comparison/ (per_run_scores.tsv, all_methods_*.csv, node_disjoint_per_seed.tsv,
ablation_auroc_crosscheck.tsv, disease_bootstrap_auroc_ci.tsv). The node-disjoint arm's metrics, CIs and agreement
are computed here with the same code (same bootstrap draws) that reproduces the stored tables for the other methods (asserted).

    python figures/main/fig6_method_comparison.py [out.png]     (B bootstrap draws: FIG6_B, default 2000)
Writes fig5_method_comparison.png/.pdf (manuscript file name), fig6b_ci_pair_vs_disease_bootstrap.tsv and
fig5_method_comparison_values.tsv next to this script.
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.transforms
from scipy.stats import rankdata

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(REPO, "src"))
from disatlas.figstyle import apply_style, panel_label, ANNOT_SIZE, LEGEND_SIZE, TICK_SIZE, DPI
D = os.path.join(REPO, "data", "derived", "fig6_comparison")
OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "fig5_method_comparison.png")

R = pd.read_csv(f"{D}/per_run_scores.tsv", sep="\t")
MT = pd.read_csv(f"{D}/all_methods_ranking_metrics.csv")
AG = pd.read_csv(f"{D}/all_methods_run_agreement.csv")
RC = pd.read_csv(f"{D}/all_methods_run_agreement_summary.csv")
AT = pd.read_csv(f"{D}/all_methods_auroc_table.csv").set_index("method")
JT = pd.read_csv(f"{D}/all_methods_top10_jaccard.csv", index_col=0)
ND_ = pd.read_csv(f"{D}/node_disjoint_per_seed.tsv", sep="\t")
assert (ND_.blind_id.values == R.blind_id.values).all() and (ND_.y.values == R.y.values).all()
R = pd.concat([R, ND_.drop(columns=["blind_id", "y"])], axis=1)

y = R.y.values.astype(int)
ND = "DisAtlas node-disjoint"
SY = ["DisAtlas", ND, "Claude-Chat Opus 5.5", "Claude-Sci Opus 5.5", "Phylo", "Gemini Pro 3.1", "TxGNN published ckpt"]
C = {"DisAtlas": "#1f4e79", ND: "#5b9bd5", "Claude-Chat Opus 5.5": "#c0504d", "Claude-Sci Opus 5.5": "#e7969c",
     "Phylo": "#2a7f62", "Gemini Pro 3.1": "#7f7f7f", "TxGNN published ckpt": "#8c564b"}
NM = {"TxGNN published ckpt": "TxGNN", "DisAtlas": "DisAtlas (full)", ND: "DisAtlas (benchmark diseases held out)", "Phylo": "Phylo Biomni",
      "Claude-Chat Opus 5.5": "Claude (chat)", "Claude-Sci Opus 5.5": "Claude (science agent)", "Gemini Pro 3.1": "Gemini Pro"}
nm = lambda s_: NM.get(s_, s_)
NMK = {"TxGNN published ckpt": "TxGNN (published ckpt)", "Phylo": "Phylo Biomni"}  # keys of the exported tables only (never drawn)
nmk = lambda s_: NMK.get(s_, s_)
wrap = lambda s_: {"TxGNN published ckpt": "TxGNN", "Claude-Sci Opus 5.5": "Claude\n(science\nagent)",
                   "DisAtlas": "DisAtlas\n(full)", ND: "DisAtlas\n(benchmark\ndiseases\nheld out)"}.get(s_, nm(s_).replace(" ", "\n", 1))
runs = {s_: [c_ for c_ in R.columns if c_.startswith(s_ + "|")] for s_ in SY}


def raw(s_):
    return R[runs[s_]].mean(axis=1).values if s_ != "DisAtlas" else R["DisAtlas"].values


def ens(s_):                       # unscored pairs ranked last
    v = raw(s_); return np.where(np.isnan(v), np.nanmin(v) - 1e-9, v)


def auc(y, s):
    r = rankdata(s); n1 = y.sum(); n0 = len(y) - n1; return (r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)


def roc(y, s):                     # tie-aware ROC
    o = np.argsort(-s, kind="mergesort"); y = y[o]; s = s[o]
    d = np.r_[np.where(np.diff(s))[0], len(s) - 1]; tp = np.cumsum(y)[d]; fp = (d + 1) - tp
    return np.r_[0, fp / fp[-1]], np.r_[0, tp / tp[-1]]


def topk_weights(v, k):            # ties at the cut-off shared fractionally
    x = np.where(np.isnan(v), -np.inf, v); w = np.zeros(len(x)); c = np.sort(x)[::-1][k - 1]
    above = x > c; tie = x == c; w[above] = 1; w[tie] = (k - above.sum()) / tie.sum(); return w


S = {s_: ens(s_) for s_ in SY}
A = {s_: auc(y, S[s_]) for s_ in SY}
# held-out (node-disjoint) model: the 7 pairs it cannot score are dropped, not ranked last 
okN = ~np.isnan(raw(ND)); yN = y[okN]; SN = raw(ND)[okN]
assert (~okN).sum() == 7, (~okN).sum()
A[ND] = auc(yN, SN)
order = sorted(SY, key=lambda s_: -A[s_])
K = int(round(0.1 * len(y)))
W_ = {s_: topk_weights(raw(s_), K) for s_ in order}
J = np.array([[np.minimum(W_[a], W_[b]).sum() / np.maximum(W_[a], W_[b]).sum() for b in order] for a in order])

# ---- value checks against the exported tables (the numbers in the current figure) ----
assert len(y) == 4364 and y.sum() == 2182 and K == 436
for s_ in [s_ for s_ in SY if s_ not in (ND, "DisAtlas")]:
    assert abs(A[s_] - MT[(MT.method == s_) & (MT.metric == "AUROC")].value.item()) < 1e-6, s_
    kA = nmk(s_) if nmk(s_) in AT.index else s_
    assert abs(A[s_] - AT.loc[kA, "AUROC"]) < 6e-5, (s_, A[s_], AT.loc[kA, "AUROC"])
ix = {s_: (nmk(s_) if nmk(s_) in JT.index else s_) for s_ in order if s_ not in (ND, "DisAtlas")}
for i, a in enumerate(order):
    for j, b in enumerate(order):
        if ND in (a, b) or "DisAtlas" in (a, b): continue
        assert abs(J[i, j] - JT.loc[ix[a], ix[b]]) < 6e-5, (a, b, J[i, j], JT.loc[ix[a], ix[b]])

# ---- node-disjoint arm: manuscript values, then metrics/agreement with the exported-table code ----
te = (R.pool == "test").values
print("AUROC full/held-out all", round(A['DisAtlas'], 4), round(A[ND], 4), "test", round(auc(y[te], S['DisAtlas'][te]), 4), round(auc(y[te & okN], raw(ND)[te & okN]), 4))
assert round(A["DisAtlas"], 3) == 0.832 and round(A[ND], 3) == 0.784, (A["DisAtlas"], A[ND])
_AG5 = pd.read_csv(f"{D}/ablation_auroc_crosscheck.tsv", sep="\t")
_a5 = _AG5[_AG5.case_set == "all"].set_index("arm").auroc
assert abs(A["DisAtlas"] - _a5["A"]) < 1e-4 and abs(A[ND] - _a5["nodedisjoint_rmref"]) < 1e-4, (A["DisAtlas"], A[ND], _a5["A"], _a5["nodedisjoint_rmref"])
print("ablation cross-check full", round(_a5["A"], 4), "held-out", round(_a5["nodedisjoint_rmref"], 4))  # full / node-disjoint, IxIDN-Matched pairs (node-disjoint: 4,357 scorable)


def exp_prec_at_k(yv, s, k):
    o = np.sort(s)[::-1]; t = o[k - 1]; gt = s > t; eq = s == t
    need = k - gt.sum(); return (yv[gt].sum() + (yv[eq].mean() * need if eq.any() else 0)) / k


def ap(yv, s):
    df_ = pd.DataFrame({"s": s, "y": yv}).groupby("s", sort=True)["y"].agg(["sum", "count"]).iloc[::-1]
    tp = df_["sum"].cumsum().values; n = df_["count"].cumsum().values; P_ = yv.sum()
    return float(np.sum((tp / n) * (df_["sum"].values / P_)))


def metrics(yv, s):
    kp = int(yv.sum()); kk = int(round(0.1 * len(yv)))
    return {"AUROC": auc(yv, s), "AUPRC": ap(yv, s), "F1@k=#positives": exp_prec_at_k(yv, s, kp),
            "Precision@top10%": exp_prec_at_k(yv, s, kk), "NPV@bottom10%": exp_prec_at_k(1 - yv, -s, kk)}


rng = np.random.default_rng(4707); BI = [rng.integers(0, len(y), len(y)) for _ in range(500)]
rows = []
for s_ in ["DisAtlas", ND]:          # DisAtlas recomputed as a check of the method
    if s_ == ND:
        rgN = np.random.default_rng(4707); m = metrics(yN, SN)
        bs = pd.DataFrame([metrics(yN[i], SN[i]) for i in (rgN.integers(0, len(yN), len(yN)) for _ in range(500))])
    else:
        m = metrics(y, S[s_]); bs = pd.DataFrame([metrics(y[i], S[s_][i]) for i in BI])
    for k_, v in m.items(): rows.append(dict(method=s_, metric=k_, value=v, lo=bs[k_].quantile(.025), hi=bs[k_].quantile(.975)))
MN = pd.DataFrame(rows)
chk = MN[MN.method == "DisAtlas"].set_index("metric")[["value", "lo", "hi"]]
ref = MT[MT.method == "DisAtlas"].set_index("metric").loc[chk.index][["value", "lo", "hi"]]
MT = pd.concat([MT[MT.method != "DisAtlas"], MN], ignore_index=True)

# ---- 95% CIs by disease-cluster bootstrap (each disease resampled
# with all its pairs, pair weight w_a*w_b, seed 4707; the figure uses B = FIG6_B draws, default 2,000; the reproduction
# check below uses the stored table's own 1,000 draws). The weight is an integer multiplicity, so each metric
# is computed on the pairs repeated by their weight; one set of draws is shared by all methods (paired).
CI_MODE = os.environ.get("FIG6_CI", "disease")
MT_pair = MT.copy()
if CI_MODE == "disease":
    nodes = pd.unique(R[["a", "b"]].values.ravel()); nix = {d_: i for i, d_ in enumerate(nodes)}
    ai_, bi_ = R.a.map(nix).values, R.b.map(nix).values
    def disease_boot(mask, B=int(os.environ.get("FIG6_B", "2000")), seed=4707):
        rg_ = np.random.default_rng(seed); out_ = []
        for _ in range(B):
            w_ = np.bincount(rg_.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); out_.append((w_[ai_] * w_[bi_])[mask])
        return out_
    # validation: stored disease-bootstrap AUROC CIs (4,357 pairs scorable by both DisAtlas models, own node set)
    _g1 = pd.read_csv(f"{D}/disease_bootstrap_auroc_ci.tsv", sep="\t")
    _g1 = _g1[_g1.scope == "all"].set_index("system")
    from sklearn.metrics import roc_auc_score
    Rv = R[okN].reset_index(drop=True); nv = pd.unique(Rv[["a", "b"]].values.ravel()); vix = {d_: i for i, d_ in enumerate(nv)}
    av, bv = Rv.a.map(vix).values, Rv.b.map(vix).values; rgv = np.random.default_rng(4707); bsv = {"full": [], "nd": []}
    for _ in range(1000):
        w_ = np.bincount(rgv.integers(0, len(nv), len(nv)), minlength=len(nv)); ww_ = (w_[av] * w_[bv]).astype(float); q_ = ww_ > 0
        bsv["full"].append(roc_auc_score(yN[q_], S["DisAtlas"][okN][q_], sample_weight=ww_[q_])); bsv["nd"].append(roc_auc_score(yN[q_], SN[q_], sample_weight=ww_[q_]))
    for k_, nm_ in [("full", "DisAtlas (full)"), ("nd", "DisAtlas (benchmark diseases held out of training, node-disjoint)")]:
        lo_, hi_ = np.percentile(bsv[k_], [2.5, 97.5]); print("CI check", k_, round(lo_, 4), round(hi_, 4), _g1.loc[nm_, ["lo", "hi"]].tolist())
        assert abs(lo_ - _g1.loc[nm_, "lo"]) < 6e-4 and abs(hi_ - _g1.loc[nm_, "hi"]) < 6e-4, (k_, lo_, hi_)
    allm = np.ones(len(y), bool); WB = disease_boot(allm)
    rows = []
    for s_ in SY:
        if s_ == ND: yv_, sv_, Wl = yN, SN, [w_[okN] for w_ in WB]
        else: yv_, sv_, Wl = y, S[s_], WB
        bs = pd.DataFrame([metrics(np.repeat(yv_, w_), np.repeat(sv_, w_)) for w_ in Wl if len(set(np.repeat(yv_, w_))) == 2])
        m = metrics(yv_, sv_)
        for k_, v in m.items(): rows.append(dict(method=s_, metric=k_, value=v, lo=bs[k_].quantile(.025), hi=bs[k_].quantile(.975), n_boot=len(bs)))
    MT = pd.DataFrame(rows)
    # value check: point estimates are unchanged
    _c = MT.merge(MT_pair, on=["method", "metric"], suffixes=("", "_p")); assert len(_c) == len(MT) and np.allclose(_c.value, _c.value_p, atol=6e-4), _c[["method", "metric", "value", "value_p"]]
    _c["width_pair"] = _c.hi_p - _c.lo_p; _c["width_disease"] = _c.hi - _c.lo
    _c[["method", "metric", "value", "lo_p", "hi_p", "lo", "hi", "width_pair", "width_disease", "n_boot"]].rename(
        columns={"lo_p": "lo_pair_bootstrap", "hi_p": "hi_pair_bootstrap", "lo": "lo_disease_bootstrap", "hi": "hi_disease_bootstrap"}).round(4).to_csv(
        os.path.join(os.path.dirname(OUT), "fig6b_ci_pair_vs_disease_bootstrap.tsv"), sep="\t", index=False)
    print("fig6b CI widths ratio disease/pair median", round(float((_c.width_disease / _c.width_pair).median()), 2))

from scipy.stats import spearmanr
def agree(s_):
    cs = runs[s_]; out_ = []
    for i in range(len(cs)):
        for j in range(i + 1, len(cs)):
            a_, b_ = R[cs[i]].values, R[cs[j]].values; ok = ~np.isnan(a_) & ~np.isnan(b_)
            out_.append(dict(method=s_, run_a=cs[i].split("|")[1], run_b=cs[j].split("|")[1], spearman=spearmanr(a_[ok], b_[ok])[0]))
    return pd.DataFrame(out_)
AG = pd.concat([AG[AG.method != "DisAtlas"], agree("DisAtlas"), agree(ND)], ignore_index=True)
print("node-disjoint", {k_: round(v, 3) for k_, v in MT[MT.method == ND].set_index("metric").value.items()},
      "CI (" + CI_MODE + " bootstrap)", {k_: (round(r_.lo, 3), round(r_.hi, 3)) for k_, r_ in MT[MT.method == ND].set_index("metric").iterrows()},
      "rho", AG[AG.method == ND].spearman.describe()[["min", "50%", "max"]].round(3).to_dict(),
      "jaccard", {nm(b): round(J[order.index(ND), j], 3) for j, b in enumerate(order)})
RC = [s_ for s_ in order if len(runs[s_]) > 1]

apply_style()
fig = plt.figure(figsize=(15.3, 12.9))
gs = fig.add_gridspec(2, 100, left=0.075, right=0.985, top=0.945, bottom=0.085, hspace=0.38, wspace=0)
a1 = fig.add_subplot(gs[0, :40]); a4 = fig.add_subplot(gs[0, 46:]); a3 = fig.add_subplot(gs[1, :40]); a2 = fig.add_subplot(gs[1, 54:])

# a: ROC
for s_ in order:
    x, t = roc(yN, SN) if s_ == ND else roc(y, S[s_]); a1.plot(x, t, color=C[s_], lw=2.2 if s_ in ("DisAtlas", ND) else 1.4, ls="--" if s_ == ND else "-", label=f"{nm(s_)}  {A[s_]:.3f}")
a1.plot([0, 1], [0, 1], ":", color="grey", lw=1); a1.set_xlim(0, 1); a1.set_ylim(0, 1.01)
a1.set_xlabel("False positive rate"); a1.set_ylabel("True positive rate")
a1.legend(title="AUROC", title_fontsize=LEGEND_SIZE, frameon=False, loc="lower right")
a1.set_title("ROC curves, IxIDN-Matched (n = 4,364; 50% positive)", loc="left")
panel_label(a1, "a")

# d: Jaccard heatmap
im = a2.imshow(J, cmap="Blues", vmin=0, vmax=1, aspect="auto")
for i_ in range(len(order)):
    for j_ in range(len(order)):
        a2.text(j_, i_, f"{J[i_, j_]:.2f}", ha="center", va="center", fontsize=TICK_SIZE,
                color="white" if J[i_, j_] > 0.55 else "black")
lab_b = [wrap(s_) for s_ in order]
a2.set_xticks(range(len(order))); a2.set_xticklabels(lab_b, fontsize=10)
a2.set_yticks(range(len(order))); a2.set_yticklabels(lab_b, fontsize=10)
a2.tick_params(length=0); [sp.set_visible(False) for sp in a2.spines.values()]
cb = fig.colorbar(im, ax=a2, fraction=0.04, pad=0.02); cb.set_label("Jaccard index"); cb.ax.tick_params(labelsize=TICK_SIZE)

# c: run-to-run agreement
rng_j = np.random.default_rng(1)
for i, s_ in enumerate(RC):
    g = AG[AG.method == s_].spearman.values
    print("fig6c run pairs", s_, len(g))
    if len(g) >= 5:   # too few run pairs for a box (Gemini Pro has 3): points only
        bp = a3.boxplot(g, positions=[i], widths=0.5, patch_artist=True, showfliers=False,
                        medianprops=dict(color="k", lw=2), whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2))
        for b in bp["boxes"]: b.set_facecolor(C[s_]); b.set_alpha(0.55); b.set_edgecolor("k")
    a3.scatter(i + rng_j.uniform(-0.12, 0.12, len(g)), g, s=26, color=C[s_], edgecolor="k", lw=0.5, zorder=3)
a3.set_xticks(range(len(RC))); a3.set_xticklabels([wrap(s_) for s_ in RC], fontsize=10)
a3.set_ylim(0, 1.12); a3.set_yticks(np.arange(0, 1.01, 0.2)); a3.set_ylabel("Spearman ρ between two runs")
a3.text(0.5, 0.99, "one point = one pair of runs; no box under 5 pairs; TxGNN (single run) not shown", transform=a3.transAxes,
        ha="center", va="top", fontsize=ANNOT_SIZE, color="0.35")
a3.set_title("Run-to-run agreement of pair scores", loc="left"); panel_label(a3, "c")
a3.grid(axis="y", color="#eee", lw=0.6); a3.set_axisbelow(True)

# b: ranking-metric small multiples inside a4's box
a4.set_title("Ranking metrics (dot = value, bar = 95% CI, diseases resampled)" if CI_MODE == "disease" else "Ranking metrics (dot = value, bar = 95% CI, pair bootstrap)", loc="left")
panel_label(a4, "b", xy=(-0.10 * a1.get_position().width / a4.get_position().width, 1.03))
a4.set_xticks([]); a4.set_yticks([]); [sp.set_visible(False) for sp in a4.spines.values()]
fig.canvas.draw()
pos = a4.get_position(); a4.set_facecolor("none")
mets = ["AUPRC", "F1@k=#positives", "Precision@top10%", "NPV@bottom10%"]
mlab = {"AUPRC": "AUPRC", "F1@k=#positives": "F1 at\nk = #pos.", "Precision@top10%": "Precision,\ntop 10%", "NPV@bottom10%": "NPV,\nbottom 10%"}
namew = 0.075; gap = 0.022; x0 = pos.x0 + namew; Wd = (pos.x1 - 0.012 - x0 - 3 * gap) / 4
yy = np.arange(len(order))[::-1]; DAX = []
p3 = a1.get_position()
for k, m in enumerate(mets):
    ax_ = fig.add_axes([x0 + k * (Wd + gap), p3.y0, Wd, p3.y1 - p3.y0 - 0.045]); DAX.append(ax_)
    g = MT[MT.metric == m].set_index("method").loc[order]
    lo_, hi_ = min(0.5, g.lo.min()), g.hi.max(); sp_ = hi_ - lo_
    for yv_, s_ in zip(yy, order):
        r_ = g.loc[s_]
        ax_.errorbar(r_.value, yv_, xerr=[[r_.value - r_.lo], [r_.hi - r_.value]], fmt="o", color=C[s_], ms=5.5,
                     mec="k", mew=0.5, elinewidth=2.0, capsize=3.5, capthick=1.4, zorder=3)
        ax_.text(r_.hi + 0.03 * sp_, yv_, f"{r_.value:.3f}", va="center", ha="left", fontsize=ANNOT_SIZE)
    ax_.axvline(0.5, color="#555", ls="--", lw=1)
    ax_.set_xlim(lo_ - 0.06 * sp_, hi_ + 0.45 * sp_); ax_.set_ylim(-0.6, len(order) - 0.4)
    tk_ = [v for v in (0.5, 0.7, 0.9) if v <= hi_ + 0.4 * sp_]; ax_.set_xticks(tk_); ax_.set_xticklabels([f"{v:g}" for v in tk_])
    ax_.set_title(mlab[m], fontsize=LEGEND_SIZE, loc="left")
    ax_.set_yticks(yy); ax_.set_yticklabels([wrap(s_) for s_ in order] if k == 0 else [])
    if k > 0: ax_.tick_params(axis="y", length=0)
    ax_.grid(axis="x", color="#eee", lw=0.6); ax_.set_axisbelow(True)
DAX[0].text(0.504, -0.45, "random", ha="left", va="bottom", fontsize=ANNOT_SIZE, color="#555")

a2.set_title("Overlap of top-10% pairs (Jaccard; 436 of 4,364 per method)", loc="left")
fig.canvas.draw()
# align panel-d title/letter with the left edge of its y tick labels
bb = a2.get_tightbbox(fig.canvas.get_renderer()).transformed(a2.transAxes.inverted())
a2.title.set_position((bb.x0, 1.0)); a2.title.set_ha("left")
panel_label(a2, "d", xy=(bb.x0 - 0.10 * a1.get_position().width / a2.get_position().width, 1.03))
fig.savefig(OUT, dpi=DPI, bbox_inches="tight"); fig.savefig(OUT[:-4] + ".pdf", bbox_inches="tight"); plt.close(fig)
print(OUT, {nm(s_): round(A[s_], 3) for s_ in order})
from disatlas.parity import write_values
_V = {f"auroc|{nm(s_)}": A[s_] for s_ in order}
for _r in MT.itertuples(): _V[f"{_r.metric}|{nm(_r.method) if _r.method in A else _r.method}"] = _r.value
write_values(OUT[:-4] + "_values.tsv", _V)
