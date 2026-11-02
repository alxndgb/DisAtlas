"""Maintainer-only: export post-threshold broad and specific channel tables for DisAtlas v3.2.

Not needed to reproduce any result (the exported files ship under data/raw/channels). Needs the internal
build outputs and the blinding salt (DISATLAS_BLIND_SALT or DISATLAS_BLIND_SALT_FILE).

usage: python scripts/maintainer/export_channels.py --arm-dir <arm graph dir> --broad-dir <shared broad dir>
          --gene-sources <tier-1 disease-gene table with 'sources'> --mechnode <mechanism-node edge table>
          --classes <mechanism-node classes tsv> --out data/raw/channels
"""
import argparse
import os
import pickle
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from disatlas import blinding as bl  # noqa: E402

ap = argparse.ArgumentParser()
for k in ("arm-dir", "broad-dir", "gene-sources", "mechnode", "classes", "out"):
    ap.add_argument("--" + k, required=True)
a = ap.parse_args()
salt = bl._salt()
ob, os_ = os.path.join(a.out, "broad"), os.path.join(a.out, "specific")
os.makedirs(ob, exist_ok=True)
os.makedirs(os_, exist_ok=True)
W = dict(index=False, compression="zstd")

# ------------------------------------------------------------------ broad channel
G = pickle.load(open(os.path.join(a.broad_dir, "G_broad.nx.pickle"), "rb"))
idm = {u: str(d.get("ID", u)) for u, d in G.nodes(data=True)}
E = sorted({tuple(sorted((idm[u], idm[v]))) for u, v in G.edges()})
pd.DataFrame(E, columns=["u", "v"]).to_parquet(os.path.join(ob, "broad_edges.parquet"), **W)
pd.DataFrame(sorted(set(idm.values())), columns=["disease_id"]).to_parquet(os.path.join(ob, "broad_nodes.parquet"), **W)

mech = os.path.join(a.broad_dir, "mech")
dg = pd.read_parquet(os.path.join(mech, "assoc_dis_gene.parquet.gzip"))
src = pd.read_parquet(a.gene_sources, columns=["disease_id", "gene", "sources"]).drop_duplicates(["disease_id", "gene"])
dg = dg.merge(src.rename(columns={"gene": "gene_symbol"}), on=["disease_id", "gene_symbol"], how="left")
srcs = dg.sources.fillna("").astype(str).str.upper().str.split(r"[;|,]")
public = srcs.map(lambda L: any(s.strip() in bl.PUBLIC_SOURCES for s in L))
dg["source_group"] = srcs.map(lambda L: "|".join(sorted({bl.source_group(s.strip()) if s.strip() else "tier_rescue"
                                                        for s in L})))
dg["target"] = dg.gene_symbol.astype(str)
dg.loc[~public, "target"] = bl.blind_series(dg.loc[~public, "gene_symbol"].astype(str), salt)
dg["target_type"] = "gene"
dg = dg[["disease_id", "target", "target_type", "evidence_tier", "tier_weight", "source_group"]]
dg.sort_values(["disease_id", "target"]).to_parquet(os.path.join(ob, "disease_gene.parquet"), **W)
dh = pd.read_parquet(os.path.join(mech, "assoc_dis_hpo.parquet.gzip"))
dh = dh[["disease_id", "hpo_id", "evidence_tier", "tier_weight"]]
dh.sort_values(["disease_id", "hpo_id"]).to_parquet(os.path.join(ob, "disease_phenotype.parquet"), **W)
print("broad", len(E), "nodes", len(set(idm.values())), "gene rows", len(dg), "blinded", int((~public).sum()),
      "hpo rows", len(dh))

# --------------------------------------------------------------- specific channel
sup = pd.read_parquet(os.path.join(a.arm_dir, "support.parquet"))
sup = sup[sup.keep].copy()
order = {m: i for i, m in enumerate(sorted(sup.module_id.unique()))}  # raw-id sort order, needed for exact float parity
mods = sorted(order)
tok = dict(zip(mods, [bl.blind_module(m, salt) for m in mods]))
assert len(set(tok.values())) == len(tok), "module token collision"
out = pd.DataFrame(dict(disease_id=sup.mondo_id.values, module_id=sup.module_id.map(tok).values,
                        module_order=sup.module_id.map(order).astype("int32").values,
                        scope=sup.scope.values, n_fam=sup.n_fam.astype("int16").values,
                        docs_max=sup.docs_max.values, families_group=sup.families.map(bl.family_groups).values))
# row order is kept as built: the projection's floating-point summation order follows it
out.to_parquet(os.path.join(os_, "disease_module_links.parquet"), **W)
print("specific links (kept)", len(out), "modules", len(mods), "blinded modules",
      int(sum(bl.is_blinded(t) for t in tok.values())))

C = pd.read_csv(a.classes, sep="\t")[["node", "mechanism_class"]]
C.sort_values(["mechanism_class", "node"]).to_csv(os.path.join(os_, "mechanism_node_classes.tsv"), sep="\t", index=False)

mn = pd.read_parquet(a.mechnode, columns=["disease_id", "node", "channel", "source", "evidence_count", "items"])
mn["items"] = mn["items"].fillna("").astype(str)
x = mn.assign(target=mn["items"].str.split(r"[;|,]")).explode("target")
x["target"] = x.target.str.strip()
x = x[x.target != ""]
r = x.source.map(bl.is_restricted).to_numpy()
x.loc[r, "target"] = bl.blind_series(x.loc[r, "target"].astype(str), salt)
x["source_group"] = x.source.map(bl.source_group)
x = x.rename(columns={"node": "mechanism_node", "channel": "target_type"})
x = x[["disease_id", "mechanism_node", "target", "target_type", "source_group", "evidence_count"]]
x = x.groupby(["disease_id", "mechanism_node", "target", "target_type", "source_group"], as_index=False).evidence_count.sum()
x.to_parquet(os.path.join(os_, "disease_target_mechanism.parquet"), **W)
print("mechanism target rows", len(x), "blinded", int(x.target.map(bl.is_blinded).sum()),
      "groups", x.source_group.value_counts().to_dict())
