"""Disease-disease GO and HPO LCA distances (the E1-E4 MechDist channels before thresholding).

Manuscript: Methods, MechDist edge definitions E1 (GO BP/MF/CC: gene-gene min LCA distance / branch
diameter, disease = symmetric best-match average over genes) and E4 (HPO: symmetric best-match average
of symptom LCA distances / diameter); appendix Figure SI_lca and Section SI_threshold_lca.

Inputs : data/00_raw/associations/disease_gene.parquet, disease_symptom.parquet,
         data/00_raw/genes/gene_go.parquet, data/00_raw/ontologies/ancestors/
Outputs: data/00_rebuilt_ontology/disease_go_distance.parquet   pairs sharing >= 1 gene
         data/00_rebuilt_ontology/disease_hpo_distance.parquet  pairs sharing >= 1 symptom (large)
Production settings: GO namespace weights BP 1.0 / MF 0.25 / CC 0.1 (summed), BP semi-regulatory
variant, GO hub-term mask at 2% (keep a gene's terms if all are masked), d_max from each matrix.
Usage: python code/00_ontology/build_disease_distances.py [--skip-hpo]
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import disease_distance as DD  # noqa: E402
import ontology as O  # noqa: E402

RAW = O.REPO_ROOT / "data" / "00_raw"
OUT = O.REPO_ROOT / "data" / "00_rebuilt_ontology"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-hpo", action="store_true", help="skip the 13.2M-pair HPO table")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    dg = pd.read_parquet(RAW / "associations" / "disease_gene.parquet")
    go = DD.go_disease_distances(DD.pairs_sharing_entity(dg, "gene_id"), dg,
                                 pd.read_parquet(RAW / "genes" / "gene_go.parquet"))
    go.to_parquet(OUT / "disease_go_distance.parquet", index=False)
    print("GO pairs", len(go), "| non-null", go[["dist_BP", "dist_MF", "dist_CC", "dist_GO"]].notna().sum().to_dict())
    if not args.skip_hpo:
        ds = pd.read_parquet(RAW / "associations" / "disease_symptom.parquet")
        hpo = DD.hpo_disease_distances(DD.pairs_sharing_entity(ds, "hpo_id"), ds)
        hpo.to_parquet(OUT / "disease_hpo_distance.parquet", index=False)
        print("HPO pairs", len(hpo), "| non-null", int(hpo.dist_HPO.notna().sum()))


if __name__ == "__main__":
    main()
