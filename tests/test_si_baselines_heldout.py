"""Supplementary Table S5, block "benchmark diseases held out of training" (4,357 pairs; test split 1,093).

table_si_baselines.py does not emit this block (commit ed1be87), but the final Supplementary Information prints it.
Every printed cell is recomputed here from data/derived/si_baselines/baselines_auroc_B2000.tsv (blocks held_out,
held_out_own_edges; one disease-cluster bootstrap, B=2000, gene-anchored model, cosine scores). The expected rows
are copied from the final SI text (lower-cased, en dashes and minus signs as '-')."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
T = pd.read_csv(ROOT / "data/derived/si_baselines/baselines_auroc_B2000.tsv", sep="\t")

EXPECTED = ['benchmark diseases held out of training (4,357 pairs; test split 1,093)', 'disatlas, embedding cosine|0.792 (0.754-0.828)|0.827 (0.766-0.876)|-', 'tf-idf graph svd, cosine; held-out diseases placed from 55 linked retained diseases|0.755 (0.712-0.794)|0.780 (0.712-0.837)|+0.037 (+0.023 to +0.052)', 'adamic-adar, same placement|0.734 (0.695-0.772)|0.740 (0.676-0.800)|+0.058 (+0.035 to +0.082)', 'weighted adamic-adar, same placement|0.724 (0.684-0.763)|0.727 (0.664-0.787)|+0.069 (+0.043 to +0.095)', 'common neighbors, same placement|0.724 (0.684-0.765)|0.726 (0.660-0.788)|+0.068 (+0.044 to +0.094)', 'tf-idf graph svd, cosine; held-out diseases keep their own graph edges|0.788 (0.747-0.827)|0.811 (0.744-0.865)|+0.004 (-0.022 to +0.028)', 'adamic-adar, own edges kept|0.770 (0.730-0.809)|0.773 (0.709-0.831)|+0.022 (-0.006 to +0.051)', 'common neighbors, own edges kept|0.760 (0.719-0.800)|0.761 (0.698-0.821)|+0.032 (+0.004 to +0.062)', 'jaccard, genes and hpo terms|0.645 (0.611-0.679)|0.657 (0.596-0.717)|+0.147 (+0.106 to +0.186)', 'graph scores with zero edges to the held-out disease|0.500 (constant)|0.500 (constant)|-']

# (printed label, block, score column)
ROWS = [
    ("disatlas, embedding cosine", "held_out", "nd_c"),
    ("tf-idf graph svd, cosine; held-out diseases placed from 55 linked retained diseases", "held_out", "tfidf_svd_place55"),
    ("adamic-adar, same placement", "held_out", "aa_place55"),
    ("weighted adamic-adar, same placement", "held_out", "waa_place55"),
    ("common neighbors, same placement", "held_out", "cn_place55"),
    ("tf-idf graph svd, cosine; held-out diseases keep their own graph edges", "held_out", "tfidf_graph_svd_cosine"),
    ("adamic-adar, own edges kept", "held_out", "adamic_adar"),
    ("common neighbors, own edges kept", "held_out", "common_neighbours"),
    ("jaccard, genes and hpo terms", "held_out", "allfeat_jaccard"),
]


def f(x):
    return f"{x:.3f}"


def s(x):
    return ("+" if x >= 0 else "-") + f"{abs(x):.3f}"


def row(label, block, score):
    a = T[(T.block == block) & (T.subset == "all") & (T.score == score)].iloc[0]
    t = T[(T.block == block) & (T.subset == "test") & (T.score == score)].iloc[0]
    cells = [label, f"{f(a.auc)} ({f(a.lo)}-{f(a.hi)})", f"{f(t.auc)} ({f(t.lo)}-{f(t.hi)})"]
    cells.append("-" if pd.isna(a.d) else f"{s(a.d)} ({s(a.dlo)} to {s(a.dhi)})")
    return "|".join(cells)


def test_heldout_block_matches_final_si():
    got = [row(*r) for r in ROWS]
    assert got == EXPECTED[1:10]


def test_heldout_header_and_zero_edge_row():
    assert EXPECTED[0] == "benchmark diseases held out of training (4,357 pairs; test split 1,093)"
    assert EXPECTED[10] == "graph scores with zero edges to the held-out disease|0.500 (constant)|0.500 (constant)|-"
    own = T[T.block == "held_out_own_edges"].set_index("score")
    ho = T[(T.block == "held_out") & (T.subset == "all")].set_index("score")
    for sc in own.index:  # the own-edges block is the all-pairs held_out row repeated
        assert abs(own.loc[sc, "auc"] - ho.loc[sc, "auc"]) < 1e-12
