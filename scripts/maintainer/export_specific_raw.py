"""Maintainer-only: export the post-filter raw inputs of the specific channel for DisAtlas v3.2.

Not needed to reproduce any result (the exported files ship under data/raw/channels/specific or Zenodo).
Needs the internal build inputs and the blinding salt (DISATLAS_BLIND_SALT or DISATLAS_BLIND_SALT_FILE).
Recomputes the internal module-support table from the internal inputs, asserts it equals the internal
arm-graph support.parquet, then writes three minimal tables:

  gene_set_counts.parquet        disease x gene set (REACTOME leaf pathways with 5..500 background genes and
                                 >= 1 disease gene, plus kept NAMED sets): n_genes, docs_max, families
  reactome_background.parquet    gene counts of the enrichment background: kind in {all, pathway, disease}
  direct_module_evidence.parquet kept direct modules (literature, pathology text, mechanism nodes) per
                                 disease and evidence family: docs (dated-document count; no ids, dates, text)

Literature and pathology modules are blinded (M_ tokens). Gene identifiers are not shipped: the
enrichment needs only the counts above, and gene-level rows for licence-restricted associations would
reveal the restricted gene through its public pathway membership.
module_order is the rank of the raw module id over all exported modules (needed for float parity).

usage: python scripts/maintainer/export_specific_raw.py --arm-dir <arm Sfull dir> --edges <edges_flagged>
          --mechnode <mechanism-node edge table> --out <dir>
"""
import argparse, os, re, sys
import numpy as np, pandas as pd
from scipy import sparse

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "src"))
from disatlas import blinding as bl  # noqa: E402

# internal build locations, passed as environment variables (never tracked):
#   DISATLAS_REACTOME_HELPERS  directory with the internal Reactome helper module (mechdist_reactome)
#   DISATLAS_ONTOLOGY_DIST_DIR internal ontology term-term distance build
#   DISATLAS_MECHDIST_DIR      internal MechDist build (disease_meta)
#   DISATLAS_RAW_DIR           internal raw-source directory
#   DISATLAS_GENE_INFO         NCBI gene_info table for Homo sapiens (csv)
sys.path.insert(0, os.environ["DISATLAS_REACTOME_HELPERS"])
import mechdist_reactome as mr  # noqa: E402

ap = argparse.ArgumentParser()
for k in ("arm-dir", "edges", "mechnode", "out"):
    ap.add_argument("--" + k, required=True)
a = ap.parse_args()
salt = bl._salt()
DATA = os.environ["DISATLAS_ONTOLOGY_DIST_DIR"]
MECH = os.environ["DISATLAS_MECHDIST_DIR"]
RAW = os.environ["DISATLAS_RAW_DIR"]
os.makedirs(a.out, exist_ok=True)
W = dict(index=False, compression="zstd")
PW_MIN, PW_MAX = 5, 500
NAMED_SRC = open(os.path.join(a.arm_dir, "slice_build.py")).read()
ns = {}; i0 = NAMED_SRC.index("NAMED = {"); exec(NAMED_SRC[i0:NAMED_SRC.index("}", i0) + 1], {}, ns)
NAMED = ns["NAMED"]
FAMS = ["AUTOAB", "CTD", "LIT", "MECHNODE", "METABASE", "OT", "PATH_ANAT", "PATH_CELL", "PATH_PHEN"]
fbit = {f: 1 << i for i, f in enumerate(FAMS)}
PUBLIC = {"OT": "open_targets", "LIT": "literature", "METABASE": "curated_restricted_a",
          "CTD": "curated_restricted_b", "AUTOAB": "literature_autoantibody",
          "PATH_ANAT": "pathology_text_anatomy", "PATH_CELL": "pathology_text_cell_type",
          "PATH_PHEN": "pathology_text_phenotype", "MECHNODE": "mechanism_node"}
def fams_public(m): return "|".join(sorted(PUBLIC[f] for f in FAMS if m & fbit[f]))
def family(src):
    if src.startswith("LIT"): return "LIT"
    if src.startswith("OT_"): return "OT"
    return {"METABASE": "METABASE", "CTD_marker_mechanism": "CTD", "AUTOAB_NLP_GENE": "AUTOAB",
            "PATH_NLP_ANATOMY": "PATH_ANAT", "PATH_NLP_CELLTYPE": "PATH_CELL", "PATH_NLP_PHENOTYPE": "PATH_PHEN"}[src]

# ---------------------------------------------------------------- admitted evidence (full slice)
U = set(pd.read_parquet(f"{MECH}/disease_meta.parquet.gzip").disease_id.astype(str))
E = pd.read_parquet(a.edges, columns=["disease_mondo", "disease_raw", "target", "target_type", "source",
                                      "n_papers", "doc_years", "drug_context"])
assert not E.drug_context.any()
E["disease"] = np.where(E.disease_mondo.isin(U), E.disease_mondo, np.where(E.disease_raw.isin(U), E.disease_raw, None))
Es = E[E.disease.notna()].copy()
islit = Es.source.str.startswith("LIT").values
cnt = np.array([(-1 if not (isinstance(s, str) and s) else sum(1 for x in s.split(";") if x)) for s in Es.doc_years.values])
docs = np.where(cnt >= 0, cnt, Es.n_papers.fillna(0).astype(float).values).astype(float)
Es["docs"] = np.where(islit, docs, 0.0)
Es = Es[~(islit & (cnt >= 0) & (cnt < 2))]
Es.loc[Es.source.str.startswith("PATH_NLP").values, "docs"] = 2.0
Es["fam"] = Es.source.map(family)

# ---------------------------------------------------------------- gene union and enrichment counts
gi = pd.read_csv(os.environ["DISATLAS_GENE_INFO"], usecols=["GeneID", "Symbol", "dbXrefs"], low_memory=False)
ens2sym = {}
for sym, x in zip(gi.Symbol, gi.dbXrefs.fillna("")):
    mm = re.search(r"Ensembl:(ENSG\d+)", x)
    if mm: ens2sym.setdefault(mm.group(1), sym)
Gn = Es[Es.target_type == "gene"]
comb = (Gn.assign(fm=Gn.fam.map(fbit).astype(np.int64)).groupby(["disease", "target"], as_index=False)
        .agg(docs=("docs", "max"), fmask=("fm", lambda s: int(np.bitwise_or.reduce(s.values)))))
comb = comb.rename(columns={"disease": "mondo_id", "target": "entity_id"})
comb["gene_symbol"] = comb.entity_id.str.replace("ENSEMBL:", "", regex=False).map(ens2sym)
Mv, idx, _ = mr.load_reactome_matrix(DATA); T = np.asarray(Mv).shape[0]; del Mv
gterms = mr.gene_pathway_rows(DATA, idx)
ens = mr.ensembl_to_entrez(RAW)
inv = {v: k for k, v in idx.items()}
cz = comb.copy(); cz["entrez"] = cz.entity_id.str.replace("ENSEMBL:", "", regex=False).map(ens)
cz = cz[cz.entrez.isin(gterms.keys())].drop_duplicates(["mondo_id", "entrez"])
genes = sorted(cz.entrez.unique()); gix = {g: i for i, g in enumerate(genes)}
dis = sorted(cz.mondo_id.unique()); di = {d: i for i, d in enumerate(dis)}
r_, c_ = zip(*[(gix[g], t) for g in genes for t in gterms[g]])
GT = sparse.csr_matrix((np.ones(len(r_), np.float32), (r_, c_)), shape=(len(genes), T))
K = np.asarray(GT.sum(0)).ravel(); ok = (K >= PW_MIN) & (K <= PW_MAX)
DG = sparse.csr_matrix((np.ones(len(cz), np.float32), (cz.mondo_id.map(di), cz.entrez.map(gix))), shape=(len(dis), len(genes)))
nn = np.asarray(DG.sum(1)).ravel()
# per (disease, pathway): gene count, max docs and OR of families over the disease's background genes
g2t = pd.DataFrame(dict(gi_=r_, t=c_)); g2t = g2t[ok[g2t.t.values]]
cz2 = cz.assign(gi_=cz.entrez.map(gix).values, di_=cz.mondo_id.map(di).values)
x = cz2[["di_", "gi_", "docs", "fmask"]].merge(g2t, on="gi_")
for i in range(len(FAMS)): x[f"b{i}"] = ((x.fmask.values >> i) & 1).astype(np.int8)
xa = x.groupby(["di_", "t"]).agg(n_genes=("gi_", "size"), docs_max=("docs", "max"),
                                 **{f"b{i}": (f"b{i}", "max") for i in range(len(FAMS))}).reset_index()
xa["fmask"] = sum(xa[f"b{i}"].values.astype(np.int64) << i for i in range(len(FAMS)))
Kmat = (DG @ GT).tocoo(); kk = pd.DataFrame(dict(di_=Kmat.row, t=Kmat.col, k=Kmat.data.astype(int)))
chk = xa.merge(kk, on=["di_", "t"], how="left"); assert (chk.k == chk.n_genes).all() and len(chk) == int(ok[kk.t].sum())
rea = pd.DataFrame(dict(disease_id=np.array(dis, dtype=object)[xa.di_.values], module_id="REACTOME:" + pd.Series([inv[t] for t in xa.t]).astype(str).values,
                        scope="reactome", n_genes=xa.n_genes.astype("int32").values, docs_max=xa.docs_max.values, fmask=xa.fmask.values))
bg = pd.concat([pd.DataFrame(dict(kind=["all"], id=["all"], n_genes=[len(genes)])),
                pd.DataFrame(dict(kind="pathway", id=["REACTOME:" + str(inv[t]) for t in np.flatnonzero(ok)], n_genes=K[ok].astype(int))),
                pd.DataFrame(dict(kind="disease", id=dis, n_genes=nn.astype(int)))], ignore_index=True)
bg["n_genes"] = bg.n_genes.astype("int32")
nm = pd.DataFrame([(k, g) for k, v in NAMED.items() for g in v.split()], columns=["module", "gene_symbol"])
y = comb.merge(nm, on="gene_symbol")
for i in range(len(FAMS)): y[f"b{i}"] = ((y.fmask.values >> i) & 1).astype(np.int8)
ya = y.groupby(["mondo_id", "module"]).agg(n_genes=("docs", "size"), docs_max=("docs", "max"),
                                          **{f"b{i}": (f"b{i}", "max") for i in range(len(FAMS))}).reset_index()
ya["fmask"] = sum(ya[f"b{i}"].values.astype(np.int64) << i for i in range(len(FAMS)))
nam = pd.DataFrame(dict(disease_id=ya.mondo_id.values, module_id=("NAMED:" + ya.module).values, scope="named",
                        n_genes=ya.n_genes.astype("int32").values, docs_max=ya.docs_max.values, fmask=ya.fmask.values))

# ---------------------------------------------------------------- direct modules
dn = Es[Es.target_type != "gene"]
pref = {"pathway": "LITPW:", "process": "LITGO:", "mechanism_entity": "LITMECH:", "anatomy": "PATH:", "cell_type": "PATH:", "phenotype": "PATH:"}
dn = pd.DataFrame(dict(disease_id=dn.disease.values, module_id=(dn.target_type.map(pref) + dn.target.astype(str)).values,
                       fam=dn.fam.values, docs=dn.docs.values)).groupby(["disease_id", "module_id", "fam"], as_index=False).docs.max()
mn = pd.read_parquet(a.mechnode); assert not mn.drug_context.any(); mn = mn[mn.disease_id.isin(U)]
mn = mn.assign(docs=[float(sum(1 for it in s.split(";") if it)) if isinstance(s, str) else 0.0 for s in mn.item_admissible_years])
mnm = mn.groupby(["disease_id", "node"], as_index=False).docs.sum()
dn = pd.concat([dn, pd.DataFrame(dict(disease_id=mnm.disease_id.values, module_id=("MECHNODE:" + mnm.node).values,
                                      fam="MECHNODE", docs=mnm.docs.values))], ignore_index=True)

# ---------------------------------------------------------------- check against the internal support table
d_ag = dn.assign(fm=dn.fam.map(fbit)).groupby(["disease_id", "module_id"]).agg(
    docs_max=("docs", "max"), fmask=("fm", lambda s: int(np.bitwise_or.reduce(s.values)))).reset_index()
d_ag["n_fam"] = [bin(v).count("1") for v in d_ag.fmask]
sup = pd.read_parquet(os.path.join(a.arm_dir, "support.parquet"))
sk = sup[sup.keep]
mine = pd.concat([rea.assign(n_fam=[bin(v).count("1") for v in rea.fmask]),
                  nam.assign(n_fam=[bin(v).count("1") for v in nam.fmask]), d_ag.assign(scope="direct")], ignore_index=True)
mine["keep"] = (mine.n_fam >= 2) | (mine.docs_max >= 2)
mk = mine[mine.keep & (mine.scope != "reactome")]
ref = sk[sk.scope != "reactome"]
m1 = ref.merge(mk, left_on=["mondo_id", "module_id", "scope"], right_on=["disease_id", "module_id", "scope"], how="outer", indicator=True)
assert (m1._merge == "both").all(), m1._merge.value_counts()
assert (m1.n_fam_x == m1.n_fam_y).all() and np.array_equal(m1.docs_max_x.values, m1.docs_max_y.values)
assert (m1.families == m1.fmask.map(lambda v: "|".join(f for f in FAMS if v & fbit[f]))).all()
rr = sup[sup.scope == "reactome"].merge(rea, left_on=["mondo_id", "module_id"], right_on=["disease_id", "module_id"], how="left")
assert rr.n_genes_y.notna().all() and (rr.n_genes_x == rr.n_genes_y).all() and np.array_equal(rr.docs_max_x.values, rr.docs_max_y.values)
assert (rr.families == rr.fmask.astype(int).map(lambda v: "|".join(f for f in FAMS if v & fbit[f]))).all()
print("internal support reproduced: direct+named kept", len(m1), "reactome rows", len(rr))

# ---------------------------------------------------------------- write (post-filter, blinded)
nam_k = nam[(nam.fmask.map(lambda v: bin(v).count("1")) >= 2) | (nam.docs_max >= 2)]
dk = d_ag[(d_ag.n_fam >= 2) | (d_ag.docs_max >= 2)][["disease_id", "module_id"]]
dn = dn.merge(dk, on=["disease_id", "module_id"])
allmods = sorted(set(rea.module_id) | set(nam_k.module_id) | set(dn.module_id))
order = {m: i for i, m in enumerate(allmods)}
tok = {m: bl.blind_module(m, salt) for m in set(dn.module_id)}
assert len(set(tok.values())) == len(tok)
gs = pd.concat([rea, nam_k], ignore_index=True)
gs = pd.DataFrame(dict(disease_id=gs.disease_id.values, module_id=gs.module_id.values,
                       module_order=gs.module_id.map(order).astype("int32").values, n_genes=gs.n_genes.values,
                       docs_max=gs.docs_max.values, families=gs.fmask.map(fams_public).values))
gs.sort_values(["disease_id", "module_order"]).to_parquet(os.path.join(a.out, "gene_set_counts.parquet"), **W)
bg.to_parquet(os.path.join(a.out, "reactome_background.parquet"), **W)
de = pd.DataFrame(dict(disease_id=dn.disease_id.values, module_id=dn.module_id.map(tok).values,
                       module_order=dn.module_id.map(order).astype("int32").values,
                       family=dn.fam.map(PUBLIC).values, docs=dn.docs.values))
de.sort_values(["disease_id", "module_order", "family"]).to_parquet(os.path.join(a.out, "direct_module_evidence.parquet"), **W)
print("gene_set_counts", len(gs), "background", len(bg), "direct evidence", len(de), "modules", len(allmods),
      "blinded", sum(bl.is_blinded(t) for t in tok.values()))
