"""Supplementary Tables SI_input_dis_gene_sources and SI_input_dis_symptom_sources: disease-gene and disease-symptom
evidence sources feeding the DisAtlas network build.

Counts only: they are taken on the integrated tier-1 association tables of the build (licensed and internal rows, not
shipped). Percentages are shares of the 13,080 candidate diseases and are recomputed here.

Input:  data/derived/input_sources/input_source_counts.json
Output: table_input_dis_gene_sources.md/.tsv and table_input_dis_symptom_sources.md/.tsv next to this script
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D = json.loads((ROOT / "data/derived/input_sources/input_source_counts.json").read_text())
N = D["n_candidate_diseases"]
BIKG = '[AZ BIKG]{custom-style="RevisionBlue"}'


def c(x):
    return f"{x:,}"


def pct(n):
    return f"{100 * n / N:.1f}%"


def b(s):
    return f"**{s}**"


def write(name, cols, sep, rows):
    pd.DataFrame(rows, columns=cols).to_csv(HERE / f"{name}.tsv", sep="\t", index=False)
    md = ["| " + " | ".join(cols) + " |", sep] + ["| " + " | ".join(r) + " |" for r in rows]
    (HERE / f"{name}.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))


PCT_COL = f'[Diseases, % of {N:,} candidates]{{custom-style="RevisionBlue"}}'
G = D["dis_gene"]; GS = G["per_source"]
GENE = [("OPENTARGETS", "OpenTargets"), ("GenCC", "GenCC"), ("OMIM Direct", "OMIM (direct)"), ("ORPHANET", "Orphanet"),
        ("BIKG_NLP_SCIENTIFIC_LITERATURE", f"{BIKG} (NLP-mined scientific literature)")]
rows = [[lab, c(GS[k]["n_assoc"]), c(GS[k]["n_diseases"]), c(GS[k]["n_genes"]), c(GS[k]["n_solo"]), pct(GS[k]["n_diseases"])]
        for k, lab in GENE]
rows.append([b("Combined (tier-1, after capping)"), b(c(G["n_rows"])), b(c(G["n_diseases"])), b(c(G["n_genes"])), "—",
             b(pct(G["n_diseases"]))])
write("table_input_dis_gene_sources",
      ["Source", "Associations", "Diseases", "Genes", "Associations exclusively from this source", PCT_COL],
      "|:---------------------------------------|----------------:|-------------:|------------:|-------------------:|----------------:|", rows)

S = D["dis_symptom"]; SS = S["per_source"]; PR = D["dis_symptom_propagation"]
SYM = [("BIKG_HPOA", f"{BIKG} (HPOA-derived)"), ("OMIM_phenotype_hpoa", "OMIM (phenotype HPOA)")]
rows = [[lab, c(SS[k]["n_assoc"]), c(SS[k]["n_diseases"]), c(SS[k]["n_symptoms"]), pct(SS[k]["n_diseases"])] for k, lab in SYM]
rows.append([b("Combined (tier-1, after capping)"), b(c(S["n_rows"])), b(c(S["n_diseases"])), b(c(S["n_symptoms"])),
             b(pct(S["n_diseases"]))])
rows.append(["MONDO-propagation enrichment (secondary augmentation, not an independent source)", c(PR["n_rows"]),
             c(PR["n_diseases"]), "—", pct(PR["n_diseases"])])
write("table_input_dis_symptom_sources", ["Source", "Associations", "Diseases", "Symptoms", PCT_COL],
      "|:---------------------------------------------------|-----------------:|--------------:|---------------:|-------------------:|", rows)
