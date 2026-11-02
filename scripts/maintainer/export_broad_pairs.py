"""Maintainer-only: export the broad-channel pair distances and the phenotype ontology slice.

Not needed to reproduce any result (the outputs ship under data/zenodo/channels_raw/broad and data/raw/ontology).
Needs the internal candidate-pair frame (one row per pair 'a__b' with dist_HPO, dist_GO, px_HPA, px_DOM) and hp.obo.

usage: python scripts/maintainer/export_broad_pairs.py --frame <candidate frame parquet> --hpo-obo <hp.obo>
          --zenodo data/zenodo/channels_raw/broad --ontology data/raw/ontology
"""
import argparse
import os

import pandas as pd

KEEP = ("id", "is_a", "relationship", "alt_id", "is_obsolete", "replaced_by", "consider", "namespace")

ap = argparse.ArgumentParser()
for k in ("frame", "hpo-obo", "zenodo", "ontology"):
    ap.add_argument("--" + k, required=True)
a = ap.parse_args()
os.makedirs(a.zenodo, exist_ok=True)
os.makedirs(a.ontology, exist_ok=True)

m = pd.read_parquet(a.frame, columns=["pair", "dist_HPO", "dist_GO", "px_HPA", "px_DOM"])
ab = m.pair.str.split("__", n=1, expand=True)
m.insert(0, "disease_b", ab[1].to_numpy())
m.insert(0, "disease_a", ab[0].to_numpy())
W = dict(index=False, compression="zstd")
m[["disease_a", "disease_b", "dist_GO", "px_HPA", "px_DOM"]].to_parquet(os.path.join(a.zenodo, "pair_gene_distances.parquet"), **W)
h = m[m.dist_HPO.notna()]
h[["disease_a", "disease_b", "dist_HPO"]].to_parquet(os.path.join(a.zenodo, "pair_hpo_distances.parquet"), **W)

# Ontology slice: only the fields the distance code reads (no header, dates or release tags).
out, keep = ["format-version: 1.2", ""], False
for line in open(a.hpo_obo):
    line = line.rstrip("\n")
    if line.startswith("["):
        keep = line == "[Term]"
        if keep:
            out += ["", "[Term]"]
    elif keep and line.split(": ", 1)[0] in KEEP:
        out.append(line.split(" ! ")[0].split(" {")[0].rstrip())
open(os.path.join(a.ontology, "hp.obo"), "w").write("\n".join(out) + "\n")
print(len(m), len(h))
