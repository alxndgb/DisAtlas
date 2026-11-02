#!/usr/bin/env python
"""SI figures: node-disjoint robustness (si_node_disjoint) and negative-set realism (si_negrealism).

Bootstrap: disease-cluster bootstrap as defined in Methods (diseases drawn with replacement, pair weight =
w[a]*w[b], pairs with an undrawn disease drop out), B = 2000, percentile 95% CI, paired deltas on shared replicates.

Inputs (data/derived/ixidn_robustness/):
  ixidn4364_per_pair_scores.tsv.gz      per-pair ensemble score p, node-disjoint ensemble score nd, splits, labels
  ixidn4364_ontology_distance.tsv       ORDO / MONDO distance per benchmark pair (ontology only)
  ixidn4364_semantic_similarity.tsv     HPO / GO Resnik similarity per pair
  random_pairs_ordo_distance.tsv.gz     ORDO distance of random disease pairs
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
B = int(os.environ.get("SI_B", 2000)); UNIT = "disease"
SYSP = "data/derived/ixidn_robustness/ixidn4364_per_pair_scores.tsv.gz"
ONTO = "data/derived/ixidn_robustness/ixidn4364_ontology_distance.tsv"
PS = "data/derived/ixidn_robustness/ixidn4364_semantic_similarity.tsv"
RND = "data/derived/ixidn_robustness/random_pairs_ordo_distance.tsv.gz"

D = pd.read_csv(f"{R0}/{SYSP}", sep="\t")
assert len(D) == 4364 and D.y.sum() == 2182
D["test"] = (D.split_a == "test") & (D.split_b == "test")
assert (D.test == (D.pool == "test")).all()
D["onc"] = D.either_is_cancer.astype(int) == 1
assert round(roc_auc_score(D.y, D.p), 3) == 0.832
D["ua"] = D.a; D["ub"] = D.b

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
assert round(roc_auc_score(N.y, N.nd), 3) == 0.784
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
ND.to_csv(f"{OUT}/si_node_disjoint_values.tsv", sep="\t", index=False, float_format="%.4f")

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
save(fig, OUT, "si_node_disjoint"); plt.close(fig)

# ===================================================================== 2. negative-set realism
od = pd.read_csv(f"{R0}/{ONTO}", sep="\t").replace(-1, np.nan)
ps = pd.read_csv(f"{R0}/{PS}", sep="\t", usecols=["blind_id", "hpo_resnik_bma", "go_resnik_genebma"])
H = D.merge(od[["blind_id", "ordo_dist", "mondo_dist"]], on="blind_id", how="left", validate="1:1").merge(ps, on="blind_id", how="left", validate="1:1")
assert len(H) == 4364
pos = (H.y == 1).values
truthy = lambda s: s.fillna(0).astype(str).str.lower().isin(["1", "1.0", "true"]).values
thr = dict(mondo=np.nanmedian(H.mondo_dist[pos]), ordo=np.nanmedian(H.ordo_dist[pos]), hpo=np.nanmedian(H.hpo_resnik_bma[pos]), go=np.nanmedian(H.go_resnik_genebma[pos]))
ax_ = dict(mondo=(H.mondo_dist <= thr["mondo"]).values, ordo=(H.ordo_dist <= thr["ordo"]).values,
           hpo=(H.hpo_resnik_bma >= thr["hpo"]).values, go=(H.go_resnik_genebma >= thr["go"]).values)
organ = truthy(H.shared_mondo_organsys); area = truthy(H.shared_ordo_area)
SUBS = [("All negatives", "neg", np.ones(len(H), bool)), ("Same organ system", "neg", organ), ("Same ORDO clinical area", "neg", area),
        ("Negatives close in MONDO", "neg", ax_["mondo"]), ("Negatives close in ORDO", "neg", ax_["ordo"]),
        ("Negatives similar by GO function", "neg", ax_["go"]), ("Negatives similar by HPO phenotype", "neg", ax_["hpo"]),
        ("Negatives close or similar on \u22652 axes", "neg", sum(a.astype(int) for a in ax_.values()) >= 2),
        ("Within stratum: organ system", "both", organ), ("Within stratum: ORDO area", "both", area)]
nneg_all = int((~pos).sum()); rows = []
for k, (nm, kind, m) in enumerate(SUBS):
    keep = (pos | (~pos & m)) if kind == "neg" else m
    for sc, sm in [("all", np.ones(len(H), bool)), ("test", H.test.values)]:
        sub = H[keep & sm].reset_index(drop=True); r = boot(sub, ["p"], seed=201 + 2 * k + (sc == "test"))["p"]
        nn = int((sub.y == 0).sum()); den = nneg_all if sc == "all" else int((~pos & sm).sum())
        rows.append(dict(subset=nm, kind=kind, scope=sc, n_pos=int(sub.y.sum()), n_neg=nn, pct_neg_kept=100 * nn / den,
                         auroc=r["auc"], lo=r["lo"], hi=r["hi"], B_used=r["B_used"]))
        print(nm, sc, rows[-1]["n_pos"], nn, round(r["auc"], 4), round(r["lo"], 4), round(r["hi"], 4), flush=True)
HN = pd.DataFrame(rows); HN["source"] = f"{SYSP}; {ONTO}; {PS}"; HN["bootstrap"] = f"{UNIT} resampling, B={B}"
HN.to_csv(f"{OUT}/si_negrealism_auroc.tsv", sep="\t", index=False, float_format="%.4f")
g1 = HN[(HN.scope == "all")].set_index("subset")
assert g1.loc["Negatives similar by HPO phenotype", "n_neg"] == 146 and g1.loc["Same organ system", "n_neg"] == 814
assert round(g1.loc["All negatives", "auroc"], 3) == 0.832 and round(g1.loc["Within stratum: organ system", "auroc"], 3) == 0.844

# panel a data (ontology distances; model-independent)
rnd = pd.read_csv(f"{R0}/{RND}", sep="\t"); rnd = rnd[rnd["set"] == "Random ORDO disease pairs"].ordo_dist.dropna().values
dneg = H.ordo_dist[~pos].dropna().values; dpos = H.ordo_dist[pos].dropna().values
A = pd.DataFrame([dict(set=s, n=len(v), q1=np.percentile(v, 25), median=np.median(v), q3=np.percentile(v, 75), mean=v.mean(), source=src)
                  for s, v, src in [("Random disease pairs", rnd, RND), ("IxIDN-Matched negatives", dneg, ONTO), ("IxIDN-Matched positives", dpos, ONTO)]])
A.to_csv(f"{OUT}/si_negrealism_ordo_distance.tsv", sep="\t", index=False, float_format="%.3f")
print(A.to_string(), flush=True)

fig = plt.figure(figsize=(7.9, 3.55))
axa = fig.add_axes([0.13, 0.12, 0.19, 0.80]); axb = fig.add_axes([0.61, 0.12, 0.27, 0.80])
bx = axa.boxplot([dpos, dneg, rnd], orientation="horizontal", widths=0.6, showfliers=False, patch_artist=True,
                 medianprops=dict(color="black", lw=1.2), whiskerprops=dict(lw=1.0), capprops=dict(lw=1.0), boxprops=dict(lw=1.0))
for p_, col in zip(bx["boxes"], ["#3B6FAD", "#D9822B", "#DDDDDD"]): p_.set_facecolor(col); p_.set_edgecolor("black")
axa.set_yticks([1, 2, 3]); axa.set_yticklabels([f"IxIDN-Matched\npositives\n(n = {len(dpos):,})", f"IxIDN-Matched\nnegatives\n(n = {len(dneg):,})", f"Random disease\npairs\n(n = {len(rnd):,})"])
axa.set_xlabel("ORDO distance (max-LCA)"); axa.set_xticks([10, 20])
ORDER = [s[0] for s in SUBS]; yy = np.arange(len(ORDER))[::-1].astype(float); off = 0.13
for i, nm in enumerate(ORDER):
    for sc, col, dy in [("all", "black", off), ("test", "#D9822B", -off)]:
        r = HN[(HN.subset == nm) & (HN.scope == sc)].iloc[0]
        axb.errorbar(r.auroc, yy[i] + dy, xerr=[[r.auroc - r.lo], [r.hi - r.auroc]], fmt="o", color=col, ms=4, lw=1.0, capsize=2)
    r = HN[(HN.subset == nm) & (HN.scope == "all")].iloc[0]
    txt = f"{r.n_neg:,} neg ({r.pct_neg_kept:.0f}%)" if r.kind == "neg" else f"{r.n_pos:,} / {r.n_neg:,} pairs"
    axb.annotate(txt, xy=(1.01, yy[i] + off), xycoords=("axes fraction", "data"), fontsize=6.5, va="center")
axb.axvline(0.5, color="#888888", ls="--", lw=0.9, zorder=0)
axb.axhline((yy[7] + yy[8]) / 2, color="#CCCCCC", lw=0.8, xmin=-0.0, xmax=1.0, clip_on=False)
axb.set_yticks(yy); axb.set_yticklabels(ORDER); axb.set_ylim(yy.min() - 0.6, yy.max() + 0.6)
xlo = np.floor((HN.lo.min() - 0.02) * 10) / 10; axb.set_xlim(min(xlo, 0.4), 0.95); axb.set_xticks(np.arange(max(0.3, min(xlo, 0.4)), 0.95, 0.1).round(1))
axb.set_xlabel("AUROC, DisAtlas (95% CI)")
h = [Line2D([], [], marker="o", color="black", ms=4, lw=1.0), Line2D([], [], marker="o", color="#D9822B", ms=4, lw=1.0)]
axb.legend(h, ["All pairs", "Test-split pairs"], loc="lower left", fontsize=7, frameon=False)
fig.text(0.015, 0.965, "a", fontweight="bold", fontsize=9); fig.text(0.375, 0.965, "b", fontweight="bold", fontsize=9)
ov = overlap_check(fig); assert not ov, ov
save(fig, OUT, "si_negrealism"); plt.close(fig)
json.dump(dict(B=B, unit=UNIT, thresholds={k: float(v) for k, v in thr.items()}, n_nd_pairs=len(N), inputs=[SYSP, ONTO, PS, RND]),
          open(f"{OUT}/si_node_disjoint_negrealism_values.json", "w"), indent=1)
print("DONE", flush=True)
