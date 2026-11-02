"""SI figure: mechanism-coverage check in dated time slices (every placeable approved expansion; no selection).
Inclusion rule (one rule, applied to all approved expansions of the 62 class drugs): the expansion disease is in the
graph in at least one slice before its first registered trial and is not already linked to an anchor disease there
(placement == 'placeable' in pair_lead_times.tsv). Failed indications under the same rule are the control.
Per-slice value = 5-seed mean stratum-matched percentile, max over anchor indications, slices <= first-trial year - 1,
slices where any anchor-target edge is already in the training graph excluded (reproduces first_cross_slice for 181/181 pairs).
Inputs: mechanism_timelines/tables/pair_lead_times.tsv, disatlas_slice_scores_5seed.tsv
Usage: si_mechanism_coverage_slices.py <src_dir> <out_dir>"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, json, numpy as np, pandas as pd, matplotlib as mpl
from matplotlib.patches import Rectangle
from disatlas import figstyle_si as fs
SRC, OUT = DERIVED / "si_mechanism_timelines", HERE; fs.apply_figure_style()
L = pd.read_csv(f"{SRC}/pair_lead_times.tsv", sep="\t"); S = pd.read_csv(f"{SRC}/disatlas_slice_scores_5seed.tsv", sep="\t")
P = L[L.placement == "placeable"].copy(); P["fty"] = pd.to_datetime(P.first_trial_date).dt.year
def traj(r):
    g = S[(S.kind == r.kind) & (S.drug == r.drug) & (S.target_id == r.target_id) & (S.year <= r.fty - 1) & S.scorable]
    e = g.groupby("year").in_G_train.any(); v = g.groupby("year").mean_pct_cm.max(); return v[~e]
T = {i: traj(r) for i, r in P.iterrows()}
for i, r in P.iterrows():  # reproduce the source crossing slice for every pair
    v = T[i]; cr = v[v >= 75].index.min() if (v >= 75).any() else np.nan
    assert (cr if pd.notna(cr) else -1) == (r.first_cross_slice if pd.notna(r.first_cross_slice) else -1), (r.drug, r.target_label)
P["cross"] = P.first_cross_slice.notna(); P["lc"] = [P.cross[i] and T[i].index.min() == P.first_cross_slice[i] for i in P.index]
N = {}
for k in ["expansion", "failed"]:
    d = P[P.kind == k]; N[k] = dict(n=len(d), cross=int(d.cross.sum()), first_scorable=int(d.lc.sum()), decoy_mean=round(float(d.decoy_cross_rate.mean()), 3), cross_pct=round(100 * d.cross.mean()))
assert (N["expansion"]["n"], N["expansion"]["cross"], N["expansion"]["first_scorable"]) == (94, 84, 55)
assert (N["failed"]["n"], N["failed"]["cross"], N["failed"]["first_scorable"]) == (87, 75, 41)
E = P[P.kind == "expansion"].copy()
SHORT = {"TNF receptor 1-associated periodic fever syndrome": "TNF receptor-associated periodic syndrome", "primary mediastinal large B-cell lymphoma": "primary mediastinal B-cell lymphoma"}
NAME = {"Complement (C5)": "Complement C5", "Complement (other)": "Complement, other targets"}
C1 = ["PD-1/PD-L1", "JAK"]; C2 = ["Complement (C5)", "Complement (other)", "TNF", "IL-6/IL-6R", "BTK", "IL-1", "IL-4R/IL-13", "IL-5/IL-5R", "IgE"]
assert set(C1 + C2) == set(E.mechanism_label)
cm = mpl.colormaps["cividis"]; nm = mpl.colors.Normalize(0, 100)
fig = mpl.pyplot.figure(figsize=(7.4, 9.6))
axs = [fig.add_axes([0.255, 0.215, 0.205, 0.745]), fig.add_axes([0.765, 0.215, 0.195, 0.745])]
nrow = max(len(E[E.mechanism_label.isin(C)]) + len(C) for C in (C1, C2))
for ax, C in zip(axs, (C1, C2)):
    y = 0; ticks, labs = [], []
    for m in C:
        ax.text(1994.2, y, NAME.get(m, m), fontsize=6, fontweight="bold", va="center", ha="right", clip_on=False); y += 1
        for i, r in E[E.mechanism_label == m].sort_values(["fty", "drug"]).iterrows():
            for yr, v in T[i].items(): ax.add_patch(Rectangle((yr, y - .42), 1, .84, fc=cm(nm(v)), ec="none"))
            if r.cross: ax.add_patch(Rectangle((r.first_cross_slice, y - .42), 1, .84, fc="none", ec="#d6604d", lw=.9, zorder=4))
            ax.plot(r.fty + .5, y, "o", ms=2.6, mfc="white", mec="k", mew=.6, zorder=5)
            ticks.append(y); labs.append(f"{r.drug} → {SHORT.get(r.target_label, r.target_label)}"); y += 1
    ax.set_yticks(ticks); ax.set_yticklabels(labs, fontsize=5); ax.tick_params(axis="y", length=0, pad=1)
    ax.set_ylim(nrow - .4, -.6); ax.set_xlim(1995, 2023); ax.set_xticks(range(1995, 2025, 5)); ax.tick_params(axis="x", labelsize=6)
    ax.spines[["top", "right", "left"]].set_visible(False); ax.set_xlabel("Time slice (training data up to end of year)", fontsize=6.5)
fig.text(.01, .985, "a", fontsize=9, fontweight="bold", va="top")
fig.text(.04, .985, f"DisAtlas percentile from the drug's first indication to each approved expansion, in every time slice before the expansion's first trial (all {N['expansion']['n']} expansions)", fontsize=7, va="top")
ca = fig.add_axes([0.80, 0.045, 0.14, 0.012]); cb = fig.colorbar(mpl.cm.ScalarMappable(norm=nm, cmap=cm), cax=ca, orientation="horizontal")
cb.set_label("Percentile among matched diseases", fontsize=6); cb.ax.tick_params(labelsize=5.5); cb.ax.axvline(75, color="#d6604d", lw=1)
b = fig.add_axes([0.255, 0.045, 0.40, 0.095]); fig.text(.01, .175, "b", fontsize=9, fontweight="bold", va="top")
fig.text(.04, .175, "Share of indications that reach the 75th percentile in at least one slice before their first trial", fontsize=7, va="top")
rows = [("Approved expansions", N["expansion"]), ("Failed indications", N["failed"])]
for j, (lab, d) in enumerate(rows):
    tot = 100 * d["cross"] / d["n"]; fs_ = 100 * d["first_scorable"] / d["n"]
    b.barh(j, fs_, color="#7a8fb8", height=.6); b.barh(j, tot - fs_, left=fs_, color="#1f3b73", height=.6)
    b.plot(100 * d["decoy_mean"], j, "|", ms=11, mew=1.6, color="#d6604d", zorder=5)
    b.text(tot + 1.5, j, f"{d['cross']}/{d['n']}", va="center", fontsize=6)
b.set_yticks([0, 1]); b.set_yticklabels([f"{l} (n = {d['n']})" for l, d in rows], fontsize=6); b.set_ylim(1.6, -.6)
b.set_xlim(0, 105); b.set_xlabel("Indications (%)", fontsize=6.5); b.tick_params(axis="x", labelsize=6); b.spines[["top", "right"]].set_visible(False)
from matplotlib.lines import Line2D; from matplotlib.patches import Patch
fig.legend(handles=[Patch(color="#1f3b73", label="reached later"), Patch(color="#7a8fb8", label="already at or above in the first slice that contains the disease"),
                    Line2D([], [], marker="|", ls="", ms=8, mew=1.6, color="#d6604d", label="matched decoy diseases, mean share"),
                    Line2D([], [], marker="o", ls="", ms=3, mfc="white", mec="k", label="first registered trial"),
                    Patch(fc="none", ec="#d6604d", label="first slice at or above 75th percentile")],
           loc="lower left", bbox_to_anchor=(.66, .065), fontsize=5.5, frameon=False, ncol=1)
N["png_bytes"] = fs.save(fig, OUT, "SI_mechanism_coverage_slices"); N["overlaps"] = fs.overlap_check(fig)
json.dump(N, open(f"{OUT}/SI_mechanism_coverage_slices_numbers.json", "w"), indent=1, default=str)
E[["mechanism_label", "drug", "target_id", "target_label", "anchors", "first_trial_date", "first_trial_nct", "first_cross_slice"]].to_csv(f"{OUT}/SI_mechanism_coverage_slices_pairs.tsv", sep="\t", index=False)
print(json.dumps(N, default=str))
