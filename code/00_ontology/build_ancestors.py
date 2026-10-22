"""Rebuild the GO (BP semi-regulatory, MF, CC) and HPO ancestor closures from the shipped OBO files.

Manuscript: Methods, "Ontological distance" (LCA distance within GO / HPO) and appendix
Figure SI_lca / Section SI_threshold_lca.

Inputs : data/00_raw/ontologies/hp.obo, data/zenodo/00_raw/ontologies/go-basic.obo
Outputs: data/00_rebuilt_ontology/ancestors/{key}_ancestors.parquet, {key}_nodes.csv, {key}_term_remap.csv
Check  : each rebuilt table equals the shipped data/00_raw/ontologies/ancestors/ table (row order free),
         and d_max (largest term-term LCA distance) is reported per key.
"""

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ontology as O  # noqa: E402

OUT = O.REPO_ROOT / "data" / "00_rebuilt_ontology" / "ancestors"


def _canon(df):
    return df.sort_values(list(df.columns)).reset_index(drop=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    built = dict(O.build_go())
    built["hpo"] = O.build_hpo()
    rows = []
    for key, (nodes, anc, remap) in built.items():
        col = "hpo_id" if key == "hpo" else "go_id"
        nodes_df = pd.DataFrame({col: nodes})
        anc.to_parquet(OUT / f"{key}_ancestors.parquet", index=False)
        nodes_df.to_csv(OUT / f"{key}_nodes.csv", index=False)
        remap.to_csv(OUT / f"{key}_term_remap.csv", index=False)
        ref = O.ANCESTORS_DIR
        same = (_canon(anc).equals(_canon(pd.read_parquet(ref / f"{key}_ancestors.parquet")))
                and _canon(nodes_df).equals(_canon(pd.read_csv(ref / f"{key}_nodes.csv")))
                and _canon(remap).equals(_canon(pd.read_csv(ref / f"{key}_term_remap.csv"))))
        d_max = int(O.distance_matrix(nodes, anc).max())
        rows.append(dict(key=key, n_terms=len(nodes), n_ancestor_rows=len(anc), d_max=d_max, equals_shipped=same))
    rep = pd.DataFrame(rows)
    print(rep.to_string(index=False))
    assert rep.equals_shipped.all(), "rebuilt ancestor tables differ from the shipped tables"
    rep.to_csv(OUT / "ancestor_build_summary.csv", index=False)


if __name__ == "__main__":
    main()
