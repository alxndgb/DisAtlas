"""Rebuild the broad MechDist graph (Habcd release) from data/raw.

Inputs (data/raw, small tables; data/zenodo/ontology_term_term_distance, term-distance slices):
  disease-gene  : gencc, orphanet, public_disease_gene (tier 1), opentargets direct (tier 1) and
                  relaxed rescue (tiers 2-4), bikg_habcd (the 8,343 released BIKG edges)
  disease-phenotype : hpo/disease_phenotype (+ _propagated, tier 9)
  specific_channels : gene2go, term_distance_params.json (data/raw/specific_channels); GO BP/MF/CC and HPO
                      term-distance slices (data/zenodo/ontology_term_term_distance)
  broad_channels    : gene_hpa_celltype_profile, gene_pfam
Steps: assemble the edge tables (lowest tier wins per pair; released BIKG edges only fill pairs not
covered by public evidence) -> four channel distances (disatlas.channels.broad_dist) -> rank-mean fusion
with CRank tie-break and TopN selection with connectivity repair (disatlas.channels.broad).

usage: python -m disatlas.graph.raw_build --out DIR
"""
import argparse
import hashlib
import json
import os

import numpy as np
import pandas as pd

from disatlas import paths
from disatlas.channels import broad, broad_dist

RAW, ZEN = paths.RAW, paths.ZENODO / "ontology_term_term_distance"
TIER_WEIGHT = {1: 1.0, 2: 0.6, 3: 0.35, 4: 0.15}
NS_KEY = {"biological_process": "go_bp_semi", "molecular_function": "go_mf", "cellular_component": "go_cc"}


def gene_edges(raw=RAW):
    """disease_id, gene_symbol, evidence_tier (+ tier_weight) in build order."""
    K = ["disease_id", "gene_symbol"]
    t1 = pd.concat([pd.read_parquet(os.path.join(raw, f)) for f in (
        "gencc/disease_gene.parquet", "orphanet/disease_gene.parquet",
        "public_disease_gene/disease_gene.parquet", "opentargets/disease_gene_direct.parquet")])[K]
    t1 = t1.drop_duplicates().assign(evidence_tier=1)
    rel = pd.read_parquet(os.path.join(raw, "opentargets/disease_gene_relaxed.parquet"))[K + ["evidence_tier"]]
    g = (pd.concat([t1, rel], ignore_index=True)
         .sort_values(K + ["evidence_tier"], kind="mergesort").drop_duplicates(K))
    bikg = pd.read_parquet(os.path.join(raw, "bikg_habcd/habcd_bikg_gene_edges.parquet"))
    bikg = bikg[K + ["evidence_tier"]]
    g = pd.concat([g, bikg], ignore_index=True).drop_duplicates(K).reset_index(drop=True)
    g["evidence_tier"] = g.evidence_tier.astype(int)
    return g.assign(tier_weight=g.evidence_tier.map(TIER_WEIGHT))


def phenotype_edges(raw=RAW):
    a = pd.read_parquet(os.path.join(raw, "hpo/disease_phenotype.parquet")).assign(evidence_tier=1)
    b = pd.read_parquet(os.path.join(raw, "hpo/disease_phenotype_propagated.parquet")).assign(evidence_tier=9)
    return pd.concat([a, b], ignore_index=True)


def _matrix(name):
    M = pd.read_parquet(ZEN / f"{name}_term_distance.parquet")
    assert list(M.index) == list(M.columns)
    return list(M.index), M.to_numpy()


def candidate_frame(raw=RAW):
    prm = json.load(open(os.path.join(raw, "specific_channels/term_distance_params.json")))
    dg = gene_edges(raw).rename(columns={"gene_symbol": "target"})
    dp = phenotype_edges(raw)
    go = pd.read_parquet(os.path.join(raw, "specific_channels/gene2go.parquet"))
    hpa = pd.read_parquet(os.path.join(raw, "broad_channels/gene_hpa_celltype_profile.parquet"))
    pfam = pd.read_parquet(os.path.join(raw, "broad_channels/gene_pfam.parquet"))
    go_mats = {ns: (*_matrix(k), prm["d_max"][k]) for ns, k in NS_KEY.items()}
    hpo_mat = (*_matrix("hpo"), prm["d_max"]["hpo"])
    m = broad_dist.candidates(dg, dp, go, hpa, pfam, go_mats, hpo_mat, prm["go_remap"])
    return m[["pair", "dist_HPO", "dist_GO"] + broad.PX]


def edge_signature(e):
    """sha256 over the sorted 'u__v' edge list."""
    k = sorted(f"{min(u, v)}__{max(u, v)}" for u, v in zip(e.u, e.v))
    return hashlib.sha256("\n".join(k).encode()).hexdigest()


def build(raw=RAW):
    m = candidate_frame(raw)
    f, quotas = broad.fuse(m.sort_values("pair").reset_index(drop=True))
    return broad.select(f, quotas), quotas, m


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    e, q, m = build()
    os.makedirs(a.out, exist_ok=True)
    m.to_parquet(os.path.join(a.out, "candidate_pairs.parquet"), index=False)
    e.to_parquet(os.path.join(a.out, "graph_edges.parquet"), index=False)
    s = {"n_edges": len(e), "n_nodes": len(set(e.u) | set(e.v)), "edge_sha256": edge_signature(e), "quotas": q}
    json.dump(s, open(os.path.join(a.out, "graph_signature.json"), "w"), indent=1)
    print(json.dumps(s))


if __name__ == "__main__":
    main()
