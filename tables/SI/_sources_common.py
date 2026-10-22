"""Shared logic for the two SI input-source tables."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import INPUTS, RAW  # noqa: E402

ASSOC = RAW / "associations"


def candidate_universe():
    """Pipeline candidate universe: MechDist-universe diseases sharing >=1 gene or HPO term
    with another disease (all evidence tiers, propagation included). 13,080 diseases."""
    dis = pd.read_parquet(RAW / "diseases" / "diseases.parquet")
    universe = set(dis.disease_id[dis.in_mechdist_universe.fillna(False).astype(bool)])
    shared = set()
    for f, key in [("disease_gene", "gene_id"), ("disease_symptom", "hpo_id")]:
        t = pd.read_parquet(ASSOC / f"{f}.parquet", columns=["disease_id", key])
        n = t.groupby(key).disease_id.nunique()
        shared |= set(t.disease_id[t[key].isin(n[n >= 2].index)])
    return universe & shared


def labels(table):
    lab = pd.read_csv(INPUTS / "source_display_labels.csv")
    return lab[lab.table == table].sort_values("display_order")


def per_source(table, key):
    t = pd.read_parquet(ASSOC / f"{table}_tier1_sources.parquet")
    nsrc = t.groupby(["disease_id", key]).source_code.transform("nunique")
    rows = []
    for _, lab in labels(table).iterrows():
        g = t[t.source_code == lab.source_code]
        rows.append(dict(label=lab.display_label, assoc=len(g), diseases=g.disease_id.nunique(),
                         items=g[key].nunique(),
                         solo_assoc=int((nsrc[t.source_code == lab.source_code] == 1).sum()),
                         disease_set=set(g.disease_id)))
    pairs = t.drop_duplicates(["disease_id", key])
    comb = dict(assoc=len(pairs), diseases=pairs.disease_id.nunique(), items=pairs[key].nunique(),
                disease_set=set(pairs.disease_id))
    return rows, comb
