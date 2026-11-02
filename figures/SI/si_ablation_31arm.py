"""SI figure: component ablation of DisAtlas (31 arms plus the gene-anchored specific-layer arm; five-seed ensembles per arm).
Input: EXPERIMENTAL_ablation_per_arm_all_readouts_FINAL.tsv. Panel a: paired change in test-split IxIDN-Matched AUROC
against the released model (disease-cluster bootstrap, B = 1,000). Panel b: share of the 159 FDA development pairs
recovered by the released model (score >= 0.5) that the arm no longer recovers (paired disease-cluster bootstrap).
The gene-anchored arm row is read from <gene_anchor_row_tsv> (written by score_gene_anchor.py / make_arm_row.py); its
mechanism-class controls were not run (no time-slice trainings), so it carries no JAK marker.
Usage: python si_ablation_gene_anchor.py <per_arm_tsv> <gene_anchor_row_tsv> <out_dir>"""
import sys, json, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
import sys as _sys; _sys.path.insert(0, str(REPO / "src"))
from disatlas.figstyle_si import apply_figure_style, overlap_check, save
apply_figure_style()
DV = REPO / "data" / "derived" / "ablation"
T = pd.read_csv(DV / "ablation_per_arm_readouts.tsv", sep="\t"); assert len(T) == 31
# Paired change in test-split AUROC and its 95% CI: disease-cluster bootstrap, B = 2000, shared replicates.
CI = pd.read_csv(DV / "ablation_paired_delta_vs_released_B2000.tsv", sep="\t").query("case_set == 'test'").set_index("arm")
for c in ["d_auroc", "d_auroc_ci_lo", "d_auroc_ci_hi", "ci_excludes_0"]:
    T[c] = T.arm.map(CI[c])
GR = pd.read_csv(DV / "gene_anchor_row.tsv", sep="\t"); assert len(GR) == 1
ga = pd.read_csv(DV / "gene_anchor_benchmark_auroc.tsv", sep="\t").query("set == 'IxIDN-Matched 4,364' and scope == 'test'").iloc[0]
GR["d_auroc"], GR["d_auroc_ci_lo"], GR["d_auroc_ci_hi"] = ga.delta, ga.delta_lo, ga.delta_hi
GR["ci_excludes_0"] = bool(ga.delta_lo > 0 or ga.delta_hi < 0)
T = pd.concat([T, GR], ignore_index=True); OUT = str(HERE)
LAB = {"gene_anchor": "Gene-anchored specific edges added", "spec_off": "Specific-edge layer removed", "broad_off": "Broad evidence edges removed",
       "unweighted": "All edge weights set to 1",
       "bikg_ge10": "Knowledge-graph module edges, stricter cut-off (10 documents)",
       "bikg_ge100": "Knowledge-graph module edges, strictest cut-off (100 documents)",
       "mech_none": "No mechanism nodes", "mech_M1M5": "Curated mechanism nodes only",
       "mech_auto277": "Automated mechanism nodes only (277)",
       "drop_REACTOME": "Reactome pathway modules removed", "drop_LITMECH": "Literature mechanism modules removed",
       "drop_LITGO": "Literature GO modules removed", "drop_LITPW": "Literature pathway modules removed",
       "drop_NAMED": "Named modules removed", "drop_lesion": "Lesion-site modules removed",
       "broad_noHPO": "Broad channel without HPO", "broad_noGO": "Broad channel without GO",
       "broad_noPfam": "Broad channel without Pfam", "broad_noHPA": "Broad channel without HPA",
       "broad_GOHPO": "Broad channel from GO and HPO only"}
FAM = {"RTK_oncogenic_angiogenic": "receptor tyrosine kinase / angiogenic", "antimetabolite_nucleotide": "antimetabolite / nucleotide",
       "JAK_STAT_interferon": "JAK–STAT / interferon", "endocrine_metabolic_other": "endocrine / metabolic",
       "checkpoint_costim": "immune checkpoint / co-stimulation", "proteasome": "proteasome", "integrin_adhesion": "integrin / adhesion",
       "cytokine_inflammatory": "inflammatory cytokine", "T_cell_immunosuppression": "T-cell immunosuppression",
       "B_cell_antibody": "B-cell / antibody", "cytokine_type2_IgE": "type-2 cytokine / IgE", "complement": "complement"}
for k, v in FAM.items(): LAB[f"drop_fam_{k}"] = v[0].upper() + v[1:]
GROUP = {"edge layer / weighting": "Edge layers and weighting", "BIKG substitution": "Knowledge-graph module threshold",
         "mechanism-node layer": "Mechanism-node layer", "single source removed": "Single evidence source removed",
         "mechanism family removed": "Mechanism family removed"}
ORDER = list(GROUP)
T["lab"] = T.arm.map(LAB); assert T.lab.notna().all(), T.arm[T.lab.isna()].tolist()
T["g"] = pd.Categorical(T.group, ORDER, ordered=True)
T = T.sort_values(["g", "frac_lost"], ascending=[True, False]).reset_index(drop=True)
assert len(T) == 32
y = []; pos = 0; gy = []
for g in ORDER:
    sub = T.index[T.group == g]
    gy.append((g, pos)); pos += 1
    for i in sub: y.append((i, pos)); pos += 1
    pos += 0.4
Y = np.zeros(len(T)); [Y.__setitem__(i, p) for i, p in y]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 7.6), sharey=True, gridspec_kw=dict(wspace=0.16))
sig = (T.d_auroc_ci_lo > 0) | (T.d_auroc_ci_hi < 0)
for ax in (a1, a2):
    for _, p in gy: ax.axhspan(p - .5, p + .5, color="none")
a1.axvline(0, color="0.4", lw=.7, ls="--")
a1.errorbar(T.d_auroc, Y, xerr=[T.d_auroc - T.d_auroc_ci_lo, T.d_auroc_ci_hi - T.d_auroc], fmt="none", ecolor="0.35", lw=.8)
a1.scatter(T.d_auroc[~sig], Y[~sig], s=16, facecolor="white", edgecolor="#0072B2", zorder=3, lw=1)
a1.scatter(T.d_auroc[sig], Y[sig], s=16, color="#0072B2", zorder=3)
a2.errorbar(T.frac_lost, Y, xerr=[T.frac_lost - T.frac_lost_lo, T.frac_lost_hi - T.frac_lost], fmt="none", ecolor="0.35", lw=.8)
jak = (T.gmpc_k_JAK < 5).values; na = T.gmpc_k_JAK.isna().values
a2.scatter(T.frac_lost[~jak & ~na], Y[~jak & ~na], s=16, color="#D55E00", zorder=3)
a2.scatter(T.frac_lost[na], Y[na], s=18, marker="s", facecolor="white", edgecolor="0.3", zorder=3, lw=1)
a2.scatter(T.frac_lost[jak], Y[jak], s=22, marker="D", facecolor="white", edgecolor="#D55E00", zorder=3, lw=1)
for i in T.index:
    a2.text(1.02, Y[i], f"{int(T.n_lost[i])}/159", va="center", ha="left", fontsize=6, transform=a2.get_yaxis_transform())
a1.set_yticks(Y); a1.set_yticklabels(T.lab, fontsize=6)
for g, p in gy:
    a1.text(-0.02, p, GROUP[g], transform=a1.get_yaxis_transform(), ha="right", va="center", fontsize=6.5, fontweight="bold")
a1.set_ylim(pos - .4, -1); a1.tick_params(axis="y", length=0)
a1.set_xlabel("Change in test-split AUROC (95% CI)"); a2.set_xlabel("Share of approved pairs no longer recovered (95% CI)")
a2.set_xlim(-0.02, 1.0); a2.set_xticks([0, .2, .4, .6, .8, 1.0]); a1.set_xticks([-.12, -.09, -.06, -.03, 0, .03]); a1.set_xlim(-.125, .03); fig.subplots_adjust(bottom=0.1)
a1.set_title("a  IxIDN-Matched, test split", loc="left", fontsize=7.5)
a2.set_title("b  FDA-approved repurposing pairs", loc="left", fontsize=7.5)
from matplotlib.lines import Line2D
h = [Line2D([], [], marker="o", ls="", color="#0072B2", ms=4, label="CI excludes zero"),
     Line2D([], [], marker="o", ls="", mfc="white", mec="#0072B2", ms=4, label="CI includes zero"),
     Line2D([], [], marker="o", ls="", color="#D55E00", ms=4, label="JAK mechanism check passes"),
     Line2D([], [], marker="D", ls="", mfc="white", mec="#D55E00", ms=4, label="JAK check fails (4 of 6 cases)"),
     Line2D([], [], marker="s", ls="", mfc="white", mec="0.3", ms=4, label="Mechanism-class checks not run")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=6.5, bbox_to_anchor=(0.55, 0.0))
ov = overlap_check(fig); b = save(fig, OUT, "si_ablation_31arm")
vals = dict(png_bytes=b, overlaps=ov, n_arms=len(T), n_sig_ixidn=int(sig.sum()), sig_arms=T.arm[sig].tolist(),
            jak_fail=T.arm[jak].tolist(), checks_not_run=T.arm[na].tolist(), ref_auroc=float(T.auroc_A.mode()[0]))
json.dump(vals, open(f"{OUT}/si_ablation_31arm_values.json", "w"), indent=1, default=str); print(vals)
