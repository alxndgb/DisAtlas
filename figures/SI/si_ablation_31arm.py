"""SI figure: component ablation of DisAtlas (31 arms, five-seed ensembles per arm), measured against the gene-anchored
model on the embedding cosine.
Input: data/derived/ablation/ablation_31arm_gene_anchored.tsv (one row per arm). Panel a: paired change in test-split
IxIDN-Matched AUROC against the gene-anchored model (disease-cluster bootstrap, B = 2,000, seed 4707). Panel b: share of
the 159 FDA development pairs that the arm no longer recovers (paired disease-cluster bootstrap, B = 2,000). Recovery is
read on the embedding cosine at a fixed all-pairs percentile, 96.595, the percentile to which the earlier 0.5 rule
corresponded (Supplementary Section S7); the per-arm counts are read from the shipped table, not recomputed here. The mechanism-class controls need
time-slice retrainings of every arm and were not run for these arms, so no JAK marker is drawn.
    python figures/SI/si_ablation_31arm.py [OUTDIR]
Writes si_ablation_31arm.png/.pdf and si_ablation_31arm_values.json (manuscript asset SI_ablation_cosine.png)."""
import sys, json, os
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from pathlib import Path
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "src"))
import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes (must precede savefig)
from disatlas.figstyle_si import apply_figure_style, overlap_check, save
apply_figure_style()
T = pd.read_csv(REPO / "data" / "derived" / "ablation" / "ablation_31arm_gene_anchored.tsv", sep="\t"); assert len(T) == 31
OUT = sys.argv[1] if len(sys.argv) > 1 else str(HERE)  # gene-anchored model is the reference: no added-arm row
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
assert len(T) == 31
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
a2.scatter(T.frac_lost, Y, s=16, color="#D55E00", zorder=3)  # no mechanism-class markers (checks not run for these arms)
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
     ]
fig.legend(handles=h, loc="lower center", ncol=2, frameon=False, fontsize=6.5, bbox_to_anchor=(0.55, 0.0))
ov = overlap_check(fig); b = save(fig, OUT, "si_ablation_31arm")
vals = dict(png_bytes=b, overlaps=ov, n_arms=len(T), n_sig_ixidn=int(sig.sum()), sig_arms=T.arm[sig].tolist(),
            jak_fail=T.arm[jak].tolist(), checks_not_run=T.arm[na].tolist(), ref_auroc=float(T.auroc_A.mode()[0]))
json.dump(vals, open(f"{OUT}/si_ablation_31arm_values.json", "w"), indent=1, default=str); print(vals)
