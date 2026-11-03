"""SI figure: broad (shared-graph) and specific (mechanism) evidence channels on FDA-approved repurposing pairs.
Development split only; the sealed holdout is not read (inputs are dev-only tables).
Inputs: data/derived/si_fda/{broad_only_fda_dev_results, broad_only_fda_dev_pair_scores_ens5, endpoints}.tsv
All three channel arms are compared on the same 5,192-disease universe (the nodes of the specific-only model; data/derived/si_fda/channel_arm_universe.tsv),
so the random-pair reference is identical across arms.
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, os, json
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from disatlas import figstyle_si as fs
SRC, OUT = DERIVED / "si_fda", HERE; fs.apply_figure_style()
R = pd.read_csv(f"{SRC}/broad_only_fda_dev_results.tsv", sep="\t")
P = pd.read_csv(f"{SRC}/broad_only_fda_dev_pair_scores_ens5.tsv", sep="\t")
E = pd.read_csv(f"{SRC}/endpoints.tsv", sep="\t")
assert not P.astype(str).apply(lambda c: c.str.contains("holdout", case=False)).any().any()
ARMS = [("full@broad_off_nodes", "Both channels (DisAtlas)", "#0072B2"),
        ("spec_off@broad_off_nodes", "Broad channel only", "#E69F00"),
        ("broad_off", "Specific channel only", "#009E73")]
N = {}
# ---- a: recovery at three thresholds
fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(2, 2, hspace=.62, wspace=.38, left=.09, right=.98, top=.93, bottom=.1)
a = fig.add_subplot(gs[0, 0]); b = fig.add_subplot(gs[0, 1]); c = fig.add_subplot(gs[1, 0]); d = fig.add_subplot(gs[1, 1])
MET = [("recovered", "Top 25%", .25), ("pct95", "Top 5%", .05), ("pct99", "Top 1%", .01)]
w = .26
for i, (arm, lab, col) in enumerate(ARMS):
    for j, (m, ml, ch) in enumerate(MET):
        r = R[(R.arm == arm) & (R.pair_set == "common") & (R.metric == m)].iloc[0]
        assert r.n_scored == 321
        x = j + (i - 1) * w
        a.bar(x, 100 * r.value, w, color=col, label=lab if j == 0 else None)
        a.errorbar(x, 100 * r.value, yerr=[[100 * (r.value - r.ci_lo)], [100 * (r.ci_hi - r.value)]], color="k", lw=.7, capsize=1.5)
        N[f"a_{arm}_{m}"] = [round(r.value, 3), round(r.ci_lo, 3), round(r.ci_hi, 3)]
for j, (m, ml, ch) in enumerate(MET):
    a.plot([j - 1.6 * w, j + 1.6 * w], [100 * ch] * 2, color="#555", ls=":", lw=.8)
a.set_xticks(range(3), [m[1] for m in MET]); a.set_ylim(0, 100); a.set_ylabel("Approved pairs recovered (%)")
a.set_xlabel("Rank of the approved pair among random disease pairs")
a.set_title("Neither channel alone matches both together", loc="left")
a.legend(frameon=False, fontsize=6, loc="upper right")
# ---- b: per-pair broad vs specific percentile
K = P[P.scorable__full].copy(); assert len(K) == 321
xb, ys = K["percentile__spec_off@broad_off_nodes"], K["percentile__broad_off"]; fr = K["percentile__full@broad_off_nodes"] >= 75
bb, ss = xb >= 75, ys >= 75
b.scatter(xb[fr], ys[fr], s=6, color="#0072B2", alpha=.6, lw=0, label="Recovered with both channels")
b.scatter(xb[~fr], ys[~fr], s=9, color="#D55E00", marker="x", lw=.7, label="Not recovered with both channels")
b.axvline(75, color="#777", lw=.6, ls="--"); b.axhline(75, color="#777", lw=.6, ls="--")
q = dict(both=int((bb & ss).sum()), broad_only=int((bb & ~ss).sum()), spec_only=int((~bb & ss).sum()), neither=int((~bb & ~ss).sum()))
for (k, xx, yy) in [("both", 98, 52), ("broad_only", 98, 3), ("spec_only", 2, 99), ("neither", 2, 3)]:
    b.text(xx, yy, f"n = {q[k]}", fontsize=6.5, ha="right" if xx > 50 else "left", va="top" if yy > 50 else "bottom", fontweight="bold", bbox=dict(fc="white", ec="none", pad=.6, alpha=.85), zorder=5)
b.set_xlim(0, 101); b.set_ylim(0, 103)
b.set_xlabel("Broad channel only (percentile)"); b.set_ylabel("Specific channel only (percentile)")
b.set_title("Each channel recovers pairs the other misses", loc="left")
b.legend(frameon=False, fontsize=6, loc="lower center", bbox_to_anchor=(.5, -.42), ncol=1, handletextpad=.3)
N["b_quadrants"] = q; N["n_scored_all_models"] = int(sum(q.values())); N["b_missed_by_both_channels_model"] = int((~fr).sum())
N["b_single_channel_recovered_but_full_missed"] = int((~fr & (bb | ss)).sum())
N["b_full_recovered_neither_single"] = int((fr & ~bb & ~ss).sum())
# ---- c: by pair class
CL = {"oncology_pan_target": "Shared cancer target", "tissue_shared": "Shared tissue", "unassigned": "Unassigned",
      "symptomatic_physiological": "Symptomatic", "convergent_effector": "Convergent effector",
      "near_duplicate": "Near-duplicate diseases", "pathogen_targeted": "Same pathogen"}
order = K.pair_class.value_counts().index.tolist()
for i, (arm, lab, col) in enumerate(ARMS):
    v = K.groupby("pair_class")[f"percentile__{arm}"].apply(lambda s: 100 * (s >= 75).mean()).reindex(order)
    c.barh(np.arange(len(order)) + (i - 1) * w, v.values, w, color=col)
    N[f"c_{arm}"] = v.round(1).to_dict()
c.set_yticks(range(len(order)), [f"{CL[o]} (n = {(K.pair_class == o).sum()})" for o in order]); c.invert_yaxis()
c.set_xlim(0, 100); c.set_xlabel("Approved pairs recovered, top 25% (%)")
c.set_title("In every pair class, the specific channel recovers\nat least as many pairs as the broad channel", loc="left")
N["c_spec_ge_broad_all_classes"] = bool(all(N["c_broad_off"][o] >= N["c_spec_off@broad_off_nodes"][o] for o in order))
# ---- d: dated (pre-trial) prioritisation
ML = [("disatlas_ens", "DisAtlas", "#0072B2"), ("shared_nbr", "Shared graph neighbors", "#777"), ("lit_shared", "Shared literature genes", "#777"),
      ("deg", "Disease degree", "#777"), ("decoy", "Matched decoy diseases", "#BBB")]
Ed = E[(E["pop"] == "dev") & (E.design == "D1")].set_index("method")
for i, (m, lab, col) in enumerate(ML):
    r = Ed.loc[m]; d.barh(i, 100 * r.top10, .6, color=col)
    if not np.isnan(r.top10_lo): d.errorbar(100 * r.top10, i, xerr=[[100 * (r.top10 - r.top10_lo)], [100 * (r.top10_hi - r.top10)]], color="k", lw=.7, capsize=1.5)
    N[f"d_{m}"] = [round(r.top10, 3), round(r.top10_lo, 3) if not np.isnan(r.top10_lo) else None, round(r.top10_hi, 3) if not np.isnan(r.top10_hi) else None]
d.set_yticks(range(len(ML)), [x[1] for x in ML]); d.invert_yaxis(); d.set_xlim(0, 100)
d.set_xlabel("Later indication in top 10% of candidates (%)")
d.set_title(f"Trained on evidence before the trial ({int(Ed.loc['disatlas_ens'].n)} pairs),\na shared-neighbor count ranks them higher than DisAtlas", loc="left")
N["d_n"] = int(Ed.loc["disatlas_ens"].n); N["d_groups"] = int(Ed.loc["disatlas_ens"].n_groups)
for ax, L in zip([a, b, c, d], "abcd"):
    ax.text(-.02, 1.07 if L in "ab" else 1.09, L, transform=ax.transAxes, fontsize=9, fontweight="bold", ha="right", va="bottom")
# claims asserted
assert N["a_spec_off@broad_off_nodes_recovered"][0] < N["a_broad_off_recovered"][0] < N["a_full@broad_off_nodes_recovered"][0]
assert q["spec_only"] == 61 and q["broad_only"] == 29
assert N["c_spec_ge_broad_all_classes"]
assert N["d_shared_nbr"][0] > N["d_disatlas_ens"][0]
N["png_bytes"] = fs.save(fig, OUT, "SI_fda_channels"); N["overlaps"] = fs.overlap_check(fig)
json.dump(N, open(f"{OUT}/SI_fda_channels_numbers.json", "w"), indent=1, default=str)
print(json.dumps(N, default=str))
