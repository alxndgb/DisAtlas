"""Maintainer-only: stage data/derived/si_rare_common/rc_annotation_counts.tsv.gz, the per-disease number of distinct
associated genes and HPO terms that define the annotation-matched null of the rare/common neighbour-coverage tables
(Supplementary Figure SI_rare_common; disatlas.rebuild.rc_coverage). Source: the authors' association loader
(02S_dis2vec/02S_mech_dist/mechdist_inputs.py load_associations, the graph-construction GDS products), as used by
revision/v32_rewrite/headline_swap/code/rc_enrichment_B.py. Counts only; no association rows are released here.
Run with the authors' environment (dis2vec_2/.venv), read-only on dis2vec_2:
    python scripts/maintainer/stage_rc_annotation_counts.py OUT.tsv.gz
"""
import os, sys
import pandas as pd
R0 = "/projects/rd/us/AI/kwfh644/dis2vec_2"
sys.path[:0] = [R0, f"{R0}/02S_dis2vec", f"{R0}/02S_dis2vec/02S_mech_dist"]
from mechdist_inputs import load_associations
HERE = os.path.dirname(os.path.abspath(__file__))
ids = pd.read_csv(os.path.join(HERE, "../../data/derived/si_rare_common/per_disease_first_other_class.tsv.gz"), sep="\t").disease_id.astype(str)
dg, dh, *_ = load_associations()
ng = dg.groupby(dg.disease_id.astype(str)).gene_symbol.nunique(); nh = dh.groupby(dh.disease_id.astype(str)).hpo_id.nunique()
T = pd.DataFrame(dict(disease_id=ids, n_genes=[int(ng.get(d, 0)) for d in ids], n_hpo=[int(nh.get(d, 0)) for d in ids]))
assert len(T) == 13019
T.to_csv(sys.argv[1], sep="\t", index=False)
print(T.describe().to_string())
