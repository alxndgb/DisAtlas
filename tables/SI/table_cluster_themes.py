"""Supplementary Table SI_cluster_themes: theme of each of the 55 disease clusters (K = 55, 13,021 diseases).

Each cluster is tested for a dominant MONDO body-system class and for an enriched mechanism or biological-process term,
and is assigned one theme: anatomy (body system alone), mixed (both), process (mechanism or process alone), neoplasm or
unresolved. The script also gives the counts, the threshold-sensitivity ranges and the Mann-Whitney contrasts that the
SI text and the Results quote.

Rule (defaults: body-system share 0.60, mechanism share 0.15, process share = mechanism share + 0.10):
  body system : the MONDO body-system class with the largest member share, excluding 'syndromic disease'; passes when
                share >= 0.60, fold >= 2 and Benjamini-Hochberg q < 0.05.
  mechanism   : any MONDO term under disease by molecular mechanism, inborn errors of metabolism or inborn error of
                immunity with share >= 0.15, fold >= 10 and q < 0.05 (deciding term: smallest p among passing terms).
  process     : the highest-ranked of the five top GO biological-process terms that is not under anatomical structure
                development (GO:0048856) or system process (GO:0003008); passes when share >= 0.25, fold >= 10, q < 0.05.
  neoplasm    : share of 'cancer or benign tumor' >= 0.40 (from the five top MONDO terms; 0 when absent).
  precedence  : neoplasm (flag, no mechanism pass, no body-system pass) > mixed > anatomy > process > unresolved.

Inputs (per-cluster annotation tables computed from the K = 55 labels, MONDO, GO, HPO and the disease-gene table):
  data/derived/si_cluster_themes/*.parquet, theme_reference.tsv (theme of each cluster, for a reproduction check)
Outputs next to this script: table_cluster_themes.{md,tsv}, table_cluster_themes_numbers.json
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
D = ROOT / "data/derived/si_cluster_themes"
rd = lambda f: pd.read_parquet(D / f"{f}.parquet")
A = rd("cluster_annotation_k55")
BL, ME = rd("bodysystem_by_cluster_long"), rd("mechanism_terms_by_cluster")
G5, M5 = rd("go_bp_top5_by_cluster"), rd("mondo_terms_top5_by_cluster")
REF = pd.read_csv(D / "theme_reference.tsv", sep="\t")
assert len(A) == 55 and A["size"].sum() == 13021
order = A.cluster.tolist()
FOLD_ORGAN, FOLD_MECH, FDR, NEO_MIN = 2.0, 10.0, 0.05, 0.40

BLo = BL[BL.category != "syndromic disease"]
organ = BLo.loc[BLo.groupby("cluster", sort=False).share.idxmax()].set_index("cluster")
G5 = G5.sort_values(["cluster", "rank"])
gocell = G5[~G5.organ_level.astype(bool)].groupby("cluster", sort=False).head(1).set_index("cluster")
neo5 = M5[M5.name == "cancer or benign tumor"].set_index("cluster").share
neo = pd.Series({c: float(neo5.get(c, 0.0)) for c in order})


def mech_pick(c, sm):
    g = ME[ME.cluster == c]
    ok = g[(g.share >= sm) & (g.fold >= FOLD_MECH) & (g.fdr < FDR)]
    r = (ok if len(ok) else g).sort_values(["p", "share"], ascending=[True, False], kind="stable").iloc[0]
    return bool(len(ok)), r


def classify(so=0.60, sm=0.15):
    sg = sm + 0.10
    rows = []
    for c in order:
        o = organ.loc[c]
        anat = bool((o.share >= so) and (o.fold >= FOLD_ORGAN) and (o.fdr < FDR))
        m1, mr = mech_pick(c, sm)
        g = gocell.loc[c] if c in gocell.index else None
        mg = bool(g is not None and (g.share >= sg) and (g.fold >= FOLD_MECH) and (g.fdr < FDR))
        isneo = neo[c] >= NEO_MIN
        th = ("neoplasm" if isneo and not m1 and not anat else "mixed" if anat and (m1 or mg) else
              "anatomy" if anat else "process" if (m1 or mg) else "unresolved")
        rows.append(dict(cluster=c, theme=th, anatomy_test=anat, organ_system=o.category, organ_share=o.share,
                         mechanism_test=m1, mechanism_term=mr["name"], mechanism_share=mr.share, mechanism_fold=mr.fold,
                         go_test=mg, go_term=None if g is None else g["name"],
                         go_share=np.nan if g is None else g.share, go_fold=np.nan if g is None else g.fold,
                         neoplasm_share=float(neo[c])))
    return pd.DataFrame(rows)


CLASSES = ["anatomy", "mixed", "process", "neoplasm", "unresolved"]
th = classify()
th.insert(1, "size", A.set_index("cluster").loc[th.cluster, "size"].values)
ref = REF.set_index("cluster").theme_class.loc[th.cluster].values
assert (th.theme.values == ref).all(), th.cluster[th.theme.values != ref].tolist()
size = th.set_index("cluster")["size"]
counts = lambda t: {k: int((t.theme == k).sum()) for k in CLASSES}
diseases = lambda t: {k: int(size[t.cluster[t.theme == k]].sum()) for k in CLASSES}


def term(r):
    if r.mechanism_test:
        return f"{r.mechanism_term} ({r.mechanism_share:.0%}, {r.mechanism_fold:.0f}×)"
    if r.go_test:
        return f"{r.go_term} ({r.go_share:.0%}, {r.go_fold:.0f}×; GO)"
    return "—"


rank = {k: i for i, k in enumerate(CLASSES)}
tt = th.assign(o=th.theme.map(rank), n=th.cluster.str[1:].astype(int)).sort_values(["o", "n"])
cols = ["Cluster", "Diseases", "Theme", "Dominant body-system class (share)", "Mechanism or process term (share, fold)"]
rows = []
for r in tt.itertuples():
    org = f"{r.organ_system} ({r.organ_share:.0%})" if r.anatomy_test else "—"
    t = r.theme + (f" ({r.neoplasm_share:.0%} cancer or benign tumor)" if r.theme == "neoplasm" else "")
    rows.append([r.cluster, f"{r.size:,}", t, org, term(r)])
pd.DataFrame(rows, columns=cols).to_csv(HERE / "table_cluster_themes.tsv", sep="\t", index=False)
md = ["| " + " | ".join(cols) + " |", "|:--|--:|:--|:--|:--|"] + ["| " + " | ".join(r) + " |" for r in rows]
(HERE / "table_cluster_themes.md").write_text("\n".join(md) + "\n", encoding="utf-8")

X = A.set_index("cluster").loc[th.cluster]
X["theme"] = th.theme.values
an, pr = X[X.theme == "anatomy"], X[X.theme == "process"]
mw = {c: dict(median_anatomy=round(float(an[c].median()), 1), median_process=round(float(pr[c].median()), 1),
              p=round(float(mannwhitneyu(an[c], pr[c]).pvalue), 3))
      for c in ["pairs_hpo_jaccard_fold", "pairs_gene_share_fold", "pairs_gobp_jaccard_fold"]}
grid = []
for so in (0.5, 0.6, 0.7):
    for sm in (0.10, 0.15, 0.20):
        cc = counts(classify(so, sm))
        grid.append(dict(body_system_dominated=cc["anatomy"] + cc["mixed"], process_enriched=cc["process"] + cc["mixed"]))
g = pd.DataFrame(grid)
cc, dd = counts(th), diseases(th)
out = dict(counts=cc, diseases=dd,
           body_system_dominated=cc["anatomy"] + cc["mixed"], process_enriched=cc["process"] + cc["mixed"],
           grid_body_system_min=int(g.body_system_dominated.min()), grid_body_system_max=int(g.body_system_dominated.max()),
           grid_process_min=int(g.process_enriched.min()), grid_process_max=int(g.process_enriched.max()),
           mann_whitney_anatomy_vs_process=mw,
           c44=dict(theme=th.set_index("cluster").theme["c44"], size=int(size["c44"])))
(HERE / "table_cluster_themes_numbers.json").write_text(json.dumps(out, indent=1) + "\n", encoding="utf-8")
print(json.dumps(dict(counts=cc, bs=out["body_system_dominated"], pe=out["process_enriched"])))
