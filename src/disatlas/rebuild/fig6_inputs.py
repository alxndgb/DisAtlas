"""Rebuild the Figure 6 / Table 1 input tables from released files.

    python -m disatlas.rebuild.fig6_inputs OUTDIR

per_run_scores.tsv     blind_id, a, b, y, pool, DisAtlas (5-seed mean cosine), DisAtlas|s<seed> (per-seed cosine), then
                       one column per comparator system|run. DisAtlas columns: data/zenodo/model_outputs/scores/
                       ixidn_matched_cosine_gene_anchored.parquet (c, c_s<seed>). Comparator columns:
                       data/raw/benchmarks/comparators/ixidn_matched_comparator_scores.tsv (pair_id = blind_id).
node_disjoint_per_seed.tsv
                       blind_id, y, DisAtlas node-disjoint|s<seed>: per seed, the mean over the folds that hold out an
                       endpoint of the pair of the re-placed endpoint's cosine (cos_rmref of
                       node_disjoint_placement_cosine_gene_anchored.parquet); NaN if a holding fold could not
                       place the endpoint (equals ndc_s<seed> of the matched cosine parquet).
Test: tests/test_rebuild_fig6_inputs.py compares both with data/derived/fig6_comparison/.
"""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SC = ROOT / "data/zenodo/model_outputs/scores"
CMP = ROOT / "data/raw/benchmarks/comparators/ixidn_matched_comparator_scores.tsv"
SEEDS = [4707, 4708, 4709, 4710, 4711]


def per_run_scores():
    S = pd.read_parquet(SC / "ixidn_matched_cosine_gene_anchored.parquet")
    out = S[["blind_id", "a", "b", "y", "pool"]].copy()
    out["DisAtlas"] = S.c
    for s in SEEDS:
        out[f"DisAtlas|s{s}"] = S[f"c_s{s}"]
    C = pd.read_csv(CMP, sep="\t").rename(columns={"pair_id": "blind_id"})
    # Phylo columns carry the stored per-run file name in the Figure 6 table
    C = C.rename(columns={k: k.replace("Phylo-Deepwalk-DDA-seed", "Phylo-Deepwalk-DDA-ixidn_expanded_BLIND_with_DDA_seed") + ".tsv"
                          for k in C.columns if k.startswith("Phylo|Phylo-Deepwalk-DDA-seed")})
    assert len(C) == len(out) == 4364 and set(C.blind_id) == set(out.blind_id)
    return out.merge(C, on="blind_id", how="left", validate="1:1")


def node_disjoint_per_seed():
    P = pd.read_parquet(SC / "node_disjoint_placement_cosine_gene_anchored.parquet")
    P = P[(P.held_a | P.held_b).to_numpy(bool)].astype({"blind_id": str, "seed": int, "cos_rmref": float})
    # NaN when any holding fold could not re-place the endpoint (the stored mean does not skip NaN; 7 pairs)
    M = P.groupby(["blind_id", "seed"]).cos_rmref.agg(lambda v: v.mean(skipna=False)).unstack("seed")
    S = pd.read_parquet(SC / "ixidn_matched_cosine_gene_anchored.parquet")[["blind_id", "y"]].astype({"blind_id": str})
    out = S.copy()
    for s in SEEDS:
        out[f"DisAtlas node-disjoint|s{s}"] = S.blind_id.map(M[s]).values
    return out


if __name__ == "__main__":
    o = Path(sys.argv[1]); o.mkdir(parents=True, exist_ok=True)
    per_run_scores().to_csv(o / "per_run_scores.tsv", sep="\t", index=False)
    node_disjoint_per_seed().to_csv(o / "node_disjoint_per_seed.tsv", sep="\t", index=False)
