#!/usr/bin/env python
"""SI figure: node-disjoint robustness (si_node_disjoint).

Bootstrap: disease-cluster bootstrap as defined in Methods (diseases drawn with replacement, pair weight =
w[a]*w[b], pairs with an undrawn disease drop out), B = 2000, percentile 95% CI, paired deltas on shared replicates.

Inputs (data/derived/ixidn_robustness/):
  ixidn4364_per_pair_scores.tsv.gz      benchmark pairs: labels, splits and pair covariates (no score columns)
Scores (data/zenodo/model_outputs/): five-seed gene-anchored embedding cosine c and node-disjoint cosine nd_c,
and the K = 55 disease clusters used as the bootstrap unit.
"""

import os, sys, json, numpy as np, pandas as pd
from sklearn.metrics import roc_auc_score
from pathlib import Path
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
import sys as _sys; _sys.path.insert(0, str(REPO / "src"))
from disatlas.figstyle_si import apply_figure_style, overlap_check, save
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R0 = str(REPO)
OUT = str(HERE)
B = int(os.environ.get("SI_B", 2000))
# resampling unit: "k55" (default; the K = 55 gene-anchored disease clusters, used by the figures) or
# "disease" (single diseases; the confidence intervals quoted in the main text). Non-default units write suffixed outputs.
UNIT = os.environ.get("SI_UNIT", "k55"); assert UNIT in ("k55", "disease")
SFX = "" if UNIT == "k55" else f"_{UNIT}"
SYSP = "data/derived/ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz"
COS = "data/zenodo/model_outputs/scores/ixidn_matched_cosine_gene_anchored.parquet"
K55 = "data/zenodo/model_outputs/clusters/disease_clusters_k55.parquet"

D = pd.read_csv(f"{R0}/{SYSP}", sep="\t")
# scores: five-seed gene-anchored embedding cosine (c) and node-disjoint cosine (nd_c) of the released model
_C = pd.read_parquet(f"{R0}/{COS}", columns=["blind_id", "c", "nd_c"])
D = D.merge(_C.rename(columns={"c": "p", "nd_c": "nd"}), on="blind_id", how="left", validate="1:1"); assert D.p.notna().all()
assert len(D) == 4364 and D.y.sum() == 2182
D["test"] = (D.split_a == "test") & (D.split_b == "test")
assert (D.test == (D.pool == "test")).all()
D["onc"] = D.either_is_cancer.astype(int) == 1
assert round(roc_auc_score(D.y, D.p), 3) == 0.834
_K = pd.read_parquet(f"{R0}/{K55}", columns=["disease_id", "cluster"]); lab = dict(zip(_K.disease_id, _K.cluster))
D["ua"] = D.a.map(lab); D["ub"] = D.b.map(lab); assert D.ua.notna().all() and D.ub.notna().all()
if UNIT == "disease": D["ua"], D["ub"] = D.a, D.b

def boot(sub, cols, seed, ref=None):
    """AUROC + percentile CI for each col; if ref, paired delta col - ref on the same replicates."""
    y = sub.y.values.astype(int); nodes = pd.unique(np.r_[sub.ua.values, sub.ub.values]); ix = {n: i for i, n in enumerate(nodes)}
    ai = np.array([ix[x] for x in sub.ua.values]); bi = np.array([ix[x] for x in sub.ub.values]); rg = np.random.default_rng(seed)
    bs = {c: [] for c in cols}
    for _ in range(B):
        w = np.bincount(rg.integers(0, len(nodes), len(nodes)), minlength=len(nodes)); ww = (w[ai] * w[bi]).astype(float); q = ww > 0
        if len(set(y[q])) < 2: continue
        for c in cols: bs[c].append(roc_auc_score(y[q], sub[c].values[q], sample_weight=ww[q]))
    out = {}
    for c in cols:
        b = np.array(bs[c]); out[c] = dict(auc=roc_auc_score(y, sub[c].values), lo=np.percentile(b, 2.5), hi=np.percentile(b, 97.5), B_used=len(b))
    if ref is not None:
        for c in cols:
            if c == ref: continue
            d = np.array(bs[c]) - np.array(bs[ref])
            out[c].update(d=out[c]["auc"] - out[ref]["auc"], dlo=np.percentile(d, 2.5), dhi=np.percentile(d, 97.5))
    return out

apply_figure_style(sizes=(8, 7, 7))
BLUE, ORANGE = "#4C72B0", "#DD8452"

# ===================================================================== 1. node-disjoint robustness
N = D[D.nd.notna()].reset_index(drop=True); assert len(N) == 4357 and N.y.sum() == 2181
assert round(roc_auc_score(N.y, N.nd), 3) == 0.792
GROUPS = [("All pairs", np.ones(len(N), bool)), ("Test split", N.test.values),
          ("Rare\u2013rare", (N.rare_common_pair == "rare-rare").values), ("Rare\u2013common", (N.rare_common_pair == "rare-common").values),
          ("Non-oncology", ~N.onc.values), ("Oncology", N.onc.values)]
rows = []
for k, (nm, m) in enumerate(GROUPS):
    sub = N[m].reset_index(drop=True); r = boot(sub, ["p", "nd"], seed=101 + k, ref="p")
    rows.append(dict(group=nm, n_pairs=len(sub), n_pos=int(sub.y.sum()), n_neg=int((sub.y == 0).sum()),
                     auroc_seen=r["p"]["auc"], seen_lo=r["p"]["lo"], seen_hi=r["p"]["hi"],
                     auroc_nd=r["nd"]["auc"], nd_lo=r["nd"]["lo"], nd_hi=r["nd"]["hi"],
                     delta=r["nd"]["d"], d_lo=r["nd"]["dlo"], d_hi=r["nd"]["dhi"], B_used=r["nd"]["B_used"],
                     ci_excludes_zero=bool(r["nd"]["dhi"] < 0 or r["nd"]["dlo"] > 0)))
    print(nm, {k_: round(v, 4) if isinstance(v, float) else v for k_, v in rows[-1].items()}, flush=True)
ND = pd.DataFrame(rows); ND["source"] = SYSP; ND["bootstrap"] = f"{UNIT} resampling, B={B}"
ND.to_csv(f"{OUT}/si_node_disjoint_values{SFX}.tsv", sep="\t", index=False, float_format="%.4f")

fig = plt.figure(figsize=(7.9, 3.75))
axa = fig.add_axes([0.20, 0.13, 0.42, 0.70]); axb = fig.add_axes([0.675, 0.13, 0.31, 0.70], sharey=axa)
yy = np.arange(len(ND))[::-1].astype(float); off = 0.15
for i, r in ND.iterrows():
    axa.errorbar(r.auroc_seen, yy[i] + off, xerr=[[r.auroc_seen - r.seen_lo], [r.seen_hi - r.auroc_seen]], fmt="o", color=BLUE, ms=5, lw=1.1, capsize=0)
    axa.errorbar(r.auroc_nd, yy[i] - off, xerr=[[r.auroc_nd - r.nd_lo], [r.nd_hi - r.auroc_nd]], fmt="o", color=ORANGE, ms=5, lw=1.1, capsize=0)
    axb.errorbar(r.delta, yy[i], xerr=[[r.delta - r.d_lo], [r.d_hi - r.delta]], fmt="o", color=ORANGE, ms=6, lw=1.1, capsize=0,
                 mfc=ORANGE if r.ci_excludes_zero else "white", mec=ORANGE, mew=1.2)
axa.set_yticks(yy); axa.set_yticklabels(ND.group); axa.set_ylim(yy.min() - 0.6, yy.max() + 0.6)
lo = min(ND.seen_lo.min(), ND.nd_lo.min()); hi = max(ND.seen_hi.max(), ND.nd_hi.max())
axa.set_xlim(np.floor((lo - 0.01) * 20) / 20, np.ceil((hi + 0.01) * 20) / 20)
axa.set_xlabel("AUROC (95% disease-cluster CI)")
axb.axvline(0, color="#444444", lw=0.9, zorder=0); axb.tick_params(axis="y", labelleft=False)
axb.set_xlabel("Node-disjoint minus seen (paired 95% CI)")
dl = min(ND.d_lo.min(), -0.01); dh = max(ND.d_hi.max(), 0.01); pad = 0.08 * (dh - dl); axb.set_xlim(dl - pad, dh + pad)
from matplotlib.ticker import MultipleLocator
def _ticks(ax, step):
    a, b = ax.get_xlim(); t = np.arange(np.ceil(a / step - 1e-9) * step, b + 1e-9, step); ax.set_xticks(np.round(t, 3)); ax.set_xlim(a, b)
_ticks(axa, 0.05); _ticks(axb, 0.05)
axa.set_title(r"$\bf{a}$  $\bf{AUROC\ by\ group\ of\ pairs}$", loc="left", fontsize=8)
axb.set_title(r"$\bf{b}$  $\bf{Paired\ change}$", loc="left", fontsize=8)
h = [Line2D([], [], marker="o", color=BLUE, ms=5, lw=1.1), Line2D([], [], marker="o", color=ORANGE, ms=5, lw=1.1)]
fig.legend(h, ["Benchmark diseases seen in training", "Node-disjoint at every step"], loc="upper center", ncol=2, bbox_to_anchor=(0.58, 1.0), fontsize=7, frameon=False)
ov = overlap_check(fig); assert not ov, ov
if UNIT == "k55": save(fig, OUT, "si_node_disjoint")
plt.close(fig)

json.dump(dict(B=B, unit=UNIT, n_nd_pairs=len(N), inputs=[SYSP, COS, K55]),
          open(f"{OUT}/si_node_disjoint_values{SFX}.json", "w"), indent=1)
print("DONE", flush=True)
