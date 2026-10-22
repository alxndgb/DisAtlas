"""Disease-disease GO and HPO distances: LCA term distance / d_max, best-match average (BMA)."""

from collections import defaultdict
from itertools import combinations

import numpy as np
import pandas as pd

import ontology as O

GO_NS_KEYS = {"biological_process": "go_bp_semi", "molecular_function": "go_mf", "cellular_component": "go_cc"}
GO_NS_WEIGHTS = {"biological_process": 1.0, "molecular_function": 0.25, "cellular_component": 0.1}
GO_HUB_FRAC = 0.02  # GO terms directly annotated to > 2% of disease genes are masked at distance time
NS_COL = {"biological_process": "dist_BP", "molecular_function": "dist_MF", "cellular_component": "dist_CC"}


def term_matrix(key):
    """(int16 distance matrix, {term: row}, d_max) from the shipped ancestor closure."""
    nodes, anc, remap = O.load_ancestors(key)
    M = O.distance_matrix(nodes, anc)
    return M, {t: i for i, t in enumerate(nodes)}, float(M.max()), remap


def entity_term_rows(df, ent_col, term_col, idx, remap):
    """{entity: sorted unique term rows} after remap; terms off the matrix axis are dropped."""
    out = defaultdict(set)
    for e, t in zip(df[ent_col], df[term_col]):
        i = idx.get(remap.get(t, t))
        if i is not None:
            out[e].add(i)
    return {e: np.array(sorted(v), dtype=np.int64) for e, v in out.items()}


def hub_terms(gene_go_ns, idx, remap, genes, frac=GO_HUB_FRAC):
    """Canonical terms directly annotated to more than ``frac`` of the on-axis disease genes."""
    sub = gene_go_ns[gene_go_ns["gene_id"].isin(set(genes))]
    canon = sub["go_id"].map(lambda t: remap.get(t, t))
    on = canon.map(lambda t: t in idx).to_numpy()
    gt = pd.DataFrame({"g": sub["gene_id"].to_numpy()[on], "t": canon.to_numpy()[on]}).drop_duplicates()
    f = gt.groupby("t")["g"].size() / gt["g"].nunique()
    return sorted(f.index[f > frac])


def masked_gene_rows(gene_go_ns, idx, remap, mask):
    """Gene term rows without ``mask`` terms; a gene annotated only to masked terms keeps them all."""
    full = entity_term_rows(gene_go_ns, "gene_id", "go_id", idx, remap)
    mrows = {idx[t] for t in mask}
    out = {}
    for g, r in full.items():
        keep = r[~np.isin(r, list(mrows))]
        out[g] = keep if len(keep) else r
    return out


def pairs_sharing_entity(df, ent_col):
    """Sorted unordered disease pairs (a < b) sharing at least one ``ent_col`` value."""
    inv = defaultdict(set)
    for d, e in zip(df["disease_id"], df[ent_col]):
        inv[e].add(d)
    pairs = set()
    for ds in inv.values():
        pairs.update(combinations(sorted(ds), 2))
    return sorted(pairs)


def _colmin(M, rows):
    return M[rows].min(axis=0)


def _bma(gg):
    """BMA of a distance matrix, dropping all-NaN rows / columns; NaN if nothing is scorable."""
    if np.isnan(gg).all():
        return np.nan
    r = np.nanmin(gg[~np.isnan(gg).all(axis=1)], axis=1)
    c = np.nanmin(gg[:, ~np.isnan(gg).all(axis=0)], axis=0)
    return 0.5 * (r.mean() + c.mean())


def gene_gene_min(M, gene_rows, genes):
    """(n_genes x n_genes) int16 min term distance between annotated genes."""
    flat = np.concatenate([gene_rows[g] for g in genes])
    starts = np.r_[0, np.cumsum([len(gene_rows[g]) for g in genes])[:-1]]
    G = np.empty((len(genes), len(genes)), dtype=np.int16)
    for k, g in enumerate(genes):
        G[k] = np.minimum.reduceat(_colmin(M, gene_rows[g])[flat], starts)
    return G


def go_disease_distances(pairs, dis_gene, gene_go):
    """DataFrame [disease_a, disease_b, dist_BP, dist_MF, dist_CC, dist_GO]; gene-gene = min term
    distance / d_max over unmasked terms, disease = BMA over genes; dist_GO = weighted sum of the
    namespaces present (NaN only if all three are NaN)."""
    dgenes = (dis_gene.drop_duplicates(["disease_id", "gene_id"])
              .groupby("disease_id")["gene_id"].apply(sorted).to_dict())
    out = pd.DataFrame({"disease_a": [a for a, _ in pairs], "disease_b": [b for _, b in pairs]})
    for ns, key in GO_NS_KEYS.items():
        M, idx, dmax, remap = term_matrix(key)
        sub = gene_go[gene_go["namespace"] == ns]
        mask = hub_terms(sub, idx, remap, dis_gene["gene_id"])
        grows = masked_gene_rows(sub, idx, remap, mask)
        genes = sorted(grows)
        gpos = {g: i for i, g in enumerate(genes)}
        G = gene_gene_min(M, grows, genes)
        del M
        vals = np.full(len(pairs), np.nan)
        for k, (a, b) in enumerate(pairs):
            ia = [gpos[g] for g in dgenes.get(a, []) if g in gpos]
            ib = [gpos[g] for g in dgenes.get(b, []) if g in gpos]
            if ia and ib:
                vals[k] = _bma(G[np.ix_(ia, ib)].astype(float) / dmax)
        out[NS_COL[ns]] = vals
    V = out[list(NS_COL.values())]
    w = np.array([GO_NS_WEIGHTS[ns] for ns in NS_COL])
    out["dist_GO"] = np.where(V.notna().any(axis=1), (V.fillna(0) * w).sum(axis=1), np.nan)
    return out


def hpo_disease_distances(pairs, dis_symptom):
    """DataFrame [disease_a, disease_b, dist_HPO]: symmetric BMA of HPO LCA distance / d_max.
    As in the production network, symptom terms are looked up on the matrix axis directly (no
    alt_id / replaced_by remap); the few obsolete or alternate HPO IDs are dropped."""
    M, idx, dmax, _ = term_matrix("hpo")
    drows = entity_term_rows(dis_symptom, "disease_id", "hpo_id", idx, {})
    dis = sorted(drows)
    colmin = {d: _colmin(M, drows[d]) for d in dis}
    del M
    vals = np.full(len(pairs), np.nan)
    for k, (a, b) in enumerate(pairs):
        if a in drows and b in drows:
            vals[k] = 0.5 * (colmin[b][drows[a]].mean() + colmin[a][drows[b]].mean()) / dmax
    return pd.DataFrame({"disease_a": [a for a, _ in pairs], "disease_b": [b for _, b in pairs], "dist_HPO": vals})
