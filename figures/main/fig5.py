"""Figure 5 (Results, cross-system comparison on IxIDN-Matched).

  a  ROC curves on all 4,364 IxIDN-Matched pairs (multi-run methods = mean over runs; unscored pairs ranked last)
  b  ranking metrics (AUPRC, F1 at k = #positives, precision top 10%, NPV bottom 10%; 95% CI, 500-draw pair bootstrap)
  c  run-to-run Spearman agreement of pair scores
  d  Jaccard overlap of top-10% pairs
Two DisAtlas arms: full and node-disjoint (diseases unseen in learning; 5 leave-diseases-out folds x 5 seeds).
LLM and TxGNN columns are the stored predictions used in the manuscript; no model is called here.

Inputs (repo-relative, the stored outputs used in the manuscript):
  data/08S_eval/ixidn4364_scores.tsv                   answer key + DisAtlas (full) ensemble and per-seed probabilities
  benchmarks/llm/data/disatlas_unseen_per_seed.tsv.gz  DisAtlas node-disjoint seeds
  benchmarks/llm/data/llm_run_scores.tsv.gz            Claude-Chat / Claude-Sci / Gemini stored runs
  benchmarks/txgnn/data/txgnn_pair_scores.tsv.gz       TxGNN published-checkpoint scores
  tables/inputs/comparator_run_scores.tsv              Phylo seeds 42-46
  benchmarks/llm/results/{ranking_metrics,run_agreement,top10_jaccard}.tsv
                                                       panels b-d for DisAtlas (both arms), Claude-Chat, Claude-Sci, Gemini
Panels b-d are computed here only for Phylo and TxGNN, which the LLM benchmark does not cover.
DisAtlas (full) = frozen ensemble (data/08S_eval p); DisAtlas (unseen in learning) = mean of the 5 node-disjoint seeds.
    python figures/main/fig5.py      -> figures/main/fig5_method_comparison.{png,pdf} and fig5_*.tsv tables
"""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.transforms
from scipy.stats import rankdata, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from disatlas_style import apply_style, panel_label, save, ANNOT_SIZE, LEGEND_SIZE, TICK_SIZE, REPO_ROOT

_rd = lambda p: pd.read_csv(os.path.join(REPO_ROOT, p), sep="\t")
# IxIDN-Matched answer key with the full-DisAtlas ensemble probability p
_SEEDS = ["p_s4707", "p_s11", "p_s202", "p_s3407", "p_s8191"]          # five training seeds of the full model
K_ = pd.read_csv(os.path.join(REPO_ROOT, "data/08S_eval/ixidn4364_scores.tsv"), sep="\t",
                 usecols=["blind_id", "a", "b", "y", "pool", "p"] + _SEEDS)
R = K_[["blind_id", "a", "b", "y", "pool", "p"] + _SEEDS].rename(
    columns={"p": "DisAtlas", **{c_: f"DisAtlas|{c_[2:]}" for c_ in _SEEDS}})
# DisAtlas node-disjoint arm: 5 seeds, each disease scored by a model that never saw it
U_ = _rd("benchmarks/llm/data/disatlas_unseen_per_seed.tsv.gz")
R = R.merge(U_.rename(columns={c_: f"DisAtlas node-disjoint|{c_[2:]}" for c_ in U_.columns if c_ != "blind_id"}), on="blind_id", how="left")
# LLM comparators: stored run outputs (no model is called)
L_ = _rd("benchmarks/llm/data/llm_run_scores.tsv.gz")
_LLM = {"claude_sci": "Claude-Sci Opus 5.5", "claude_chat": "Claude-Chat Opus 5.5", "gemini_pro_3_1": "Gemini Pro 3.1"}
R = R.merge(L_.rename(columns={c_: f"{_LLM[c_.rsplit('_run', 1)[0]]}|run{c_.rsplit('_run', 1)[1]}" for c_ in L_.columns if c_ != "blind_id"}),
            on="blind_id", how="left")
# TxGNN published checkpoint: stored pair scores; unmapped or collapsed pairs are unscored
T_ = _rd("benchmarks/txgnn/data/txgnn_pair_scores.tsv.gz")
T_["TxGNN published ckpt|run1"] = T_.score.where(T_.mapped & ~T_.collapsed)
R = R.merge(T_[["blind_id", "TxGNN published ckpt|run1"]], on="blind_id", how="left")
# Phylo: stored scores for five seeds
P_ = _rd("tables/inputs/comparator_run_scores.tsv")
R = R.merge(P_[["blind_id"] + [c_ for c_ in P_.columns if c_.startswith("Phylo|")]], on="blind_id", how="left")
assert len(R) == 4364 and R.blind_id.is_unique
y = R.y.values.astype(int)
ND = "DisAtlas node-disjoint"
SY = ["DisAtlas", ND, "Claude-Chat Opus 5.5", "Claude-Sci Opus 5.5", "Phylo", "Gemini Pro 3.1", "TxGNN published ckpt"]
C = {"DisAtlas": "#1f4e79", ND: "#5b9bd5", "Claude-Chat Opus 5.5": "#c0504d", "Claude-Sci Opus 5.5": "#e7969c",
     "Phylo": "#2a7f62", "Gemini Pro 3.1": "#7f7f7f", "TxGNN published ckpt": "#8c564b"}
NM = {"TxGNN published ckpt": "TxGNN (published ckpt)", "DisAtlas": "DisAtlas (full)", ND: "DisAtlas (unseen in learning)"}
nm = lambda s_: NM.get(s_, s_)
wrap = lambda s_: {"TxGNN published ckpt": "TxGNN\n(published ckpt)",
                   "DisAtlas": "DisAtlas\n(full)", ND: "DisAtlas\n(unseen in\nlearning)"}.get(s_, nm(s_).replace(" ", "\n", 1))
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
order = sorted(SY, key=lambda s_: -A[s_])
K = int(round(0.1 * len(y)))
W_ = {s_: topk_weights(raw(s_), K) for s_ in order}
J = np.array([[np.minimum(W_[a], W_[b]).sum() / np.maximum(W_[a], W_[b]).sum() for b in order] for a in order])


# ---- manuscript values ----
assert len(y) == 4364 and y.sum() == 2182 and K == 436
te = (R.pool == "test").values
assert round(A[ND], 3) == 0.742 and round(auc(y[te], S[ND][te]), 3) == 0.754, (A[ND], auc(y[te], S[ND][te]))


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


# ---- panels b-d: shipped LLM-benchmark results for the methods it covers; Phylo and TxGNN computed here ----
BM = {"DisAtlas (full)": "DisAtlas", "DisAtlas (unseen in learning)": ND, "Claude Chat (Opus 5.5)": "Claude-Chat Opus 5.5",
      "Claude Sci (Opus 5.5)": "Claude-Sci Opus 5.5", "Gemini Pro 3.1": "Gemini Pro 3.1"}
NEW = [s_ for s_ in SY if s_ not in BM.values()]                       # Phylo, TxGNN published ckpt
_bm = lambda f: _rd(f"benchmarks/llm/results/{f}.tsv")
MT_b = _bm("ranking_metrics").assign(method=lambda d: d.method.map(BM))
for s_ in BM.values():                                                   # stored AUROC equals the one recomputed from scores
    assert abs(MT_b.query("method == @s_ and metric == 'AUROC'").value.item() - A[s_]) < 5e-7, s_
rng = np.random.default_rng(4707); BI = [rng.integers(0, len(y), len(y)) for _ in range(500)]
rows = []
for s_ in NEW:
    m = metrics(y, S[s_]); bs = pd.DataFrame([metrics(y[i], S[s_][i]) for i in BI])
    for k_, v in m.items(): rows.append(dict(method=s_, metric=k_, value=v, lo=bs[k_].quantile(.025), hi=bs[k_].quantile(.975)))
MT = pd.concat([MT_b, pd.DataFrame(rows)], ignore_index=True)
MT = MT.assign(_o=MT.method.map(SY.index)).sort_values("_o", kind="stable").drop(columns="_o").reset_index(drop=True)


def agree(s_):
    cs = runs[s_]; out_ = []
    for i in range(len(cs)):
        for j in range(i + 1, len(cs)):
            a_, b_ = R[cs[i]].values, R[cs[j]].values; ok = ~np.isnan(a_) & ~np.isnan(b_)
            out_.append(dict(method=s_, run_a=cs[i].split("|")[1], run_b=cs[j].split("|")[1], spearman=spearmanr(a_[ok], b_[ok])[0]))
    return pd.DataFrame(out_)
AG = pd.concat([_bm("run_agreement").assign(method=lambda d: d.method.map(BM))]
               + [agree(s_) for s_ in NEW if len(runs[s_]) > 1], ignore_index=True)
JB = _bm("top10_jaccard").set_index("method").rename(index=BM, columns=BM)
for a_ in BM.values():
    for b_ in BM.values(): J[order.index(a_), order.index(b_)] = JB.loc[a_, b_]
MT.assign(method=MT.method.map(nm)).to_csv(os.path.join(HERE, "fig5b_ranking_metrics.tsv"), sep="\t", index=False, float_format="%.6f")
AG.assign(method=AG.method.map(nm)).to_csv(os.path.join(HERE, "fig5c_run_agreement.tsv"), sep="\t", index=False, float_format="%.6f")
pd.DataFrame(J, index=[nm(s_) for s_ in order], columns=[nm(s_) for s_ in order]).to_csv(
    os.path.join(HERE, "fig5d_top10_jaccard.tsv"), sep="\t", float_format="%.6f")
print({nm(s_): round(A[s_], 3) for s_ in order})
RC = [s_ for s_ in order if len(runs[s_]) > 1]

apply_style()
fig = plt.figure(figsize=(15.3, 12.9))
gs = fig.add_gridspec(2, 100, left=0.075, right=0.985, top=0.945, bottom=0.085, hspace=0.38, wspace=0)
a1 = fig.add_subplot(gs[0, :40]); a4 = fig.add_subplot(gs[0, 46:]); a3 = fig.add_subplot(gs[1, :40]); a2 = fig.add_subplot(gs[1, 54:])

# a: ROC
for s_ in order:
    x, t = roc(y, S[s_]); a1.plot(x, t, color=C[s_], lw=2.2 if s_ in ("DisAtlas", ND) else 1.4, ls="--" if s_ == ND else "-", label=f"{nm(s_)}  {A[s_]:.3f}")
a1.plot([0, 1], [0, 1], ":", color="grey", lw=1); a1.set_xlim(0, 1); a1.set_ylim(0, 1.01)
a1.set_xlabel("False positive rate"); a1.set_ylabel("True positive rate")
a1.legend(title="AUROC (all 4,364 pairs)", title_fontsize=LEGEND_SIZE, frameon=False, loc="lower right")
a1.set_title("ROC curves, IxIDN-Matched (n = 4,364; 50% positive)", loc="left")
panel_label(a1, "a")

# d: Jaccard heatmap
im = a2.imshow(J, cmap="Blues", vmin=0, vmax=1, aspect="auto")
for i_ in range(len(order)):
    for j_ in range(len(order)):
        a2.text(j_, i_, f"{J[i_, j_]:.2f}", ha="center", va="center", fontsize=TICK_SIZE,
                color="white" if J[i_, j_] > 0.55 else "black")
lab_b = [wrap(s_).replace("(published ckpt)", "(published\nckpt)").replace("(unseen in learning)", "(unseen in\nlearning)") for s_ in order]
a2.set_xticks(range(len(order))); a2.set_xticklabels(lab_b, fontsize=10)
a2.set_yticks(range(len(order))); a2.set_yticklabels(lab_b, fontsize=10)
a2.tick_params(length=0); [sp.set_visible(False) for sp in a2.spines.values()]
cb = fig.colorbar(im, ax=a2, fraction=0.04, pad=0.02); cb.set_label("Jaccard index"); cb.ax.tick_params(labelsize=TICK_SIZE)

# c: run-to-run agreement
rng_j = np.random.default_rng(1)
for i, s_ in enumerate(RC):
    g = AG[AG.method == s_].spearman.values
    bp = a3.boxplot(g, positions=[i], widths=0.5, patch_artist=True, showfliers=False,
                    medianprops=dict(color="k", lw=2), whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2))
    for b in bp["boxes"]: b.set_facecolor(C[s_]); b.set_alpha(0.55); b.set_edgecolor("k")
    a3.scatter(i + rng_j.uniform(-0.12, 0.12, len(g)), g, s=26, color=C[s_], edgecolor="k", lw=0.5, zorder=3)
a3.set_xticks(range(len(RC))); a3.set_xticklabels([wrap(s_) for s_ in RC], fontsize=10)
a3.set_ylim(0, 1.12); a3.set_yticks(np.arange(0, 1.01, 0.2)); a3.set_ylabel("Spearman ρ between two runs")
a3.text(0.5, 0.99, "one point = one pair of runs; TxGNN (single run) not shown", transform=a3.transAxes,
        ha="center", va="top", fontsize=ANNOT_SIZE, color="0.35")
a3.set_title("Run-to-run agreement of pair scores", loc="left"); panel_label(a3, "c")
a3.grid(axis="y", color="#eee", lw=0.6); a3.set_axisbelow(True)

# b: ranking-metric small multiples inside a4's box
a4.set_title("Ranking metrics (dot = value, bar = 95% CI, pair bootstrap)", loc="left")
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
save(fig, "fig5_method_comparison", bbox_inches="tight"); plt.close(fig)
