"""Supplementary Table SI_baselines: DisAtlas embedding cosine against training-free graph and annotation scores on the
IxIDN-Matched benchmark, with the benchmark diseases in the graph (4,364 pairs). The held-out-of-training block
(4,357 pairs) is not emitted (author decision 2026-11-10).

AUROCs and paired differences (DisAtlas minus score; all pairs and test split) come from one disease-cluster bootstrap with
2,000 draws on the gene-anchored model; MONDO SimGIC is a point estimate.

Input:  data/derived/si_baselines/baselines_auroc_B2000.tsv
Output: table_si_baselines.tsv and table_si_baselines.md next to this script
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T = pd.read_csv(ROOT / "data/derived/si_baselines/baselines_auroc_B2000.tsv", sep="\t")

# (label, block of the AUROC rows, score column, block of the difference row or None).
# Within each block (a row with block None starts a block) rows are emitted in descending order of the
# unrounded all-pairs AUROC. The whole block of benchmark diseases held out of training is not emitted (author
# decision 2026-11-10): its header, the DisAtlas held-out row, the placement and own-edge rows, the Jaccard
# genes-and-HPO row and the zero-edge row. Its values stay in baselines_auroc_B2000.tsv (blocks held_out,
# held_out_own_edges). With one block left, its header row (benchmark diseases in the graph, 4,364 pairs; test split
# 1,094) is not emitted either; the rows are the in-graph rows.
ROWS = [
    ("DisAtlas, embedding cosine", "in_graph", "c", None),
    ("TF-IDF graph SVD, cosine", "in_graph", "tfidf_graph_svd_cosine__full_graph", "in_graph"),
    ("Adjacency SVD (unweighted), cosine", "in_graph", "adjacency_svd_binary_cosine", "in_graph"),
    ("Adjacency SVD (weighted), cosine", "in_graph", "adjacency_svd_weighted_cosine", "in_graph"),
    ("Adamic–Adar", "in_graph", "adamic_adar__full_graph", "in_graph"),
    ("Common neighbors", "in_graph", "common_neighbours__full_graph", "in_graph"),
    ("Weighted Adamic–Adar", "in_graph", "weighted_adamic_adar__full_graph", "in_graph"),
    ("Degree product", "in_graph", "degree_product__full_graph", "in_graph"),
    ("Jaccard, genes and HPO terms", "in_graph", "allfeat_jaccard", "in_graph"),
    ("Jaccard, genes", "in_graph", "gene_jaccard", "in_graph"),
    ("Jaccard, HPO terms", "in_graph", "hpo_jaccard", "in_graph"),
    ("MONDO SimGIC", "in_graph", "mondo_simgic", None),
]


def g(x):
    return f"{x:.3f}".replace("-", "−")


def s(x):
    return f"{x:+.3f}".replace("-", "−")


def get(block, subset, score):
    r = T[(T.block == block) & (T.subset == subset) & (T.score == score)]
    assert len(r) == 1, (block, subset, score, len(r))
    return r.iloc[0]


def auc(r):
    return g(r.auc) if pd.isna(r.lo) else f"{g(r.auc)} ({g(r.lo)}–{g(r.hi)})"


def row(lab, blk, sc, dblk):
    a, t = get(blk, "all", sc), get(blk, "test", sc)
    dd = lambda sub: "–" if dblk is None else (lambda r: f"{s(r.d)} ({s(r.dlo)} to {s(r.dhi)})")(get(dblk, sub, sc))
    return a.auc, [lab, auc(a), auc(t), dd("all"), dd("test")]


out, block = [], []
for lab, blk, sc, dblk in ROWS + [(None, None, None, None)]:
    if blk is None:
        out += [r for _, r in sorted(block, key=lambda x: -x[0])]  # stable: ties keep listed order
        block = []
        if lab is not None:
            out.append([lab, "", "", "", ""])
        continue
    block.append(row(lab, blk, sc, dblk))
assert not any("held" in r[0].lower() for r in out)
H = ["Score", "All pairs, AUROC (95% CI)", "Test split, AUROC (95% CI)", "Difference, all pairs (95% CI)",
     "Difference, test split (95% CI)"]
pd.DataFrame(out, columns=H).to_csv(HERE / "table_si_baselines.tsv", sep="\t", index=False)
md = ["| " + " | ".join(H) + " |", "|---|---|---|---|---|"] + ["| " + " | ".join(r) + " |" for r in out]
(HERE / "table_si_baselines.md").write_text("\n".join(md) + "\n")
print(len(out))
