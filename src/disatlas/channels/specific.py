"""Specific channel: rebuild the disease-module links from the post-filter raw evidence.

Inputs (data/raw/channels/specific, or data/zenodo/channels_raw/specific for the large table):
  gene_set_counts.parquet         disease x gene set: REACTOME leaf pathways (5..500 background genes, >= 1
                                  disease gene) and NAMED sets; n_genes, docs_max, families
  reactome_background.parquet     background gene counts: kind 'all' (N), 'pathway' (K), 'disease' (n)
  direct_module_evidence.parquet  literature, pathology-text and mechanism-node modules per disease and family

Recipe: one-sided hypergeometric test per disease and pathway (k >= 2), Benjamini-Hochberg over the
pathways with k >= 1, q < 0.05 -> enriched REACTOME module; plus NAMED and direct modules. A link is kept
when it has >= 2 evidence families or >= 2 documents. graph/build.py later drops modules linked to
fewer than 2 diseases and weights a link by n_fam + ln(1 + docs).

usage: python -m disatlas.channels.specific --out DIR [--data data] [--cold]
(--cold is accepted for pipeline symmetry: this channel has no frozen part, so both modes rebuild it.)
"""
import argparse
import os

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

Q_MAX, K_MIN = 0.05, 2
GROUP = {"open_targets": "open_targets", "literature": "literature", "curated_restricted_a": "curated_restricted",
         "curated_restricted_b": "curated_restricted", "literature_autoantibody": "literature_autoantibody",
         "pathology_text_anatomy": "pathology_text", "pathology_text_cell_type": "pathology_text",
         "pathology_text_phenotype": "pathology_text", "mechanism_node": "mechanism_node"}


def raw_path(data, name):
    for d in ("raw/channels/specific", "zenodo/channels_raw/specific"):
        p = os.path.join(data, d, name)
        if os.path.exists(p):
            return p
    raise FileNotFoundError(name)


def enriched(gs, bg):
    """(disease_id, module_id) of REACTOME pathways enriched at BH q < Q_MAX."""
    r = gs[gs.module_id.str.startswith("REACTOME:")]
    N = int(bg.loc[bg.kind == "all", "n_genes"].iloc[0])
    K = r.module_id.map(bg[bg.kind == "pathway"].set_index("id").n_genes).to_numpy(float)
    n = r.disease_id.map(bg[bg.kind == "disease"].set_index("id").n_genes).to_numpy(float)
    k = r.n_genes.to_numpy(float)
    p = np.full(len(r), np.nan)
    t = k >= K_MIN  # one survival-function call per distinct (k, K, n): same values, far fewer calls
    u, inv = np.unique(np.c_[k[t], K[t], n[t]], axis=0, return_inverse=True)
    p[t] = hypergeom.sf(u[:, 0] - 1, N, u[:, 1], u[:, 2])[inv.ravel()]
    keep = np.zeros(len(r), bool)
    for _, ix in r.groupby("disease_id", sort=False).indices.items():
        tested = ix[k[ix] >= K_MIN]
        if not len(tested):
            continue
        pv = p[tested]
        o = np.argsort(pv)
        q = pv[o] * len(ix) / (np.arange(len(pv)) + 1)  # BH over all pathways with k >= 1
        q = np.minimum.accumulate(q[::-1])[::-1]
        qq = np.empty_like(q)
        qq[o] = np.minimum(q, 1)
        keep[tested[qq < Q_MAX]] = True
    return r[keep]


def build(data="data"):
    gs = pd.read_parquet(raw_path(data, "gene_set_counts.parquet"))
    bg = pd.read_parquet(raw_path(data, "reactome_background.parquet"))
    de = pd.read_parquet(raw_path(data, "direct_module_evidence.parquet"))
    rea = enriched(gs, bg).assign(scope="reactome")
    nam = gs[gs.module_id.str.startswith("NAMED:")].assign(scope="named")
    d = de.groupby(["disease_id", "module_id", "module_order"], as_index=False).agg(
        docs_max=("docs", "max"), families=("family", lambda s: "|".join(sorted(set(s)))))
    cols = ["disease_id", "module_id", "module_order", "scope", "docs_max", "families"]
    L = pd.concat([rea[cols], nam[cols], d.assign(scope="direct")[cols]], ignore_index=True)
    fam = L.families.str.split("|")
    L["n_fam"] = fam.map(len).astype("int16")
    L = L[(L.n_fam >= 2) | (L.docs_max >= 2)]
    L["families_group"] = fam[L.index].map(lambda f: "|".join(sorted({GROUP[x] for x in f})))
    # rows in (disease, raw module id) order; module_order re-ranked over the kept modules
    L = L.sort_values(["disease_id", "module_order"], kind="mergesort").reset_index(drop=True)
    L["module_order"] = L.module_order.rank(method="dense").astype("int32") - 1
    return L[["disease_id", "module_id", "module_order", "scope", "n_fam", "docs_max", "families_group"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--data", default="data")
    ap.add_argument("--cold", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    L = build(a.data)
    L.to_parquet(os.path.join(a.out, "disease_module_links.parquet"), index=False, compression="zstd")
    print("disease_module_links", len(L), "diseases", L.disease_id.nunique(), "modules", L.module_id.nunique())


if __name__ == "__main__":
    main()
