#!/usr/bin/env python
"""Figure 1 (rare and common diseases: genetic evidence and shared genes) from Open Targets.

Subcommands
-----------
tables   Build the source tables (CSV) for Figure 1a, Figure 1b and the SI per-source panel
         from an Open Targets Platform release and the frozen rare/common label inputs.
figure   Render Figure 1 (panels a-d) and the SI per-source panel from the tables.

Inputs (all paths are arguments; nothing is hard-coded)
-------------------------------------------------------
--ot-root       Open Targets Platform release directory with the bulk parquet datasets
                ``disease/``, ``association_by_datatype_direct/`` and
                ``association_by_datasource_direct/`` (see scripts/fetch_open_targets.py), or, with
                ``--input-format csv2024``, a directory with the 2024 CSV exports
                ``opentargets-diseases.csv``, ``ot_assoc_lookup.csv`` and
                ``ot_gene_association_by_datasource.csv``.
--labels-root   Directory with the frozen label inputs:
                ``mondo_rare_common.parquet`` (MONDO rare/common table),
                ``orphanet_idout.csv`` (Orphanet NON RARE / OBSOLETE list),
                ``relabel_round1.xlsx`` and ``relabel_round2.xlsx`` (expert relabel lists),
                and optionally ``id_crosswalk.csv`` (retired Open Targets id -> current id),
                ``mondo_rare_common_venn.parquet`` (MONDO table used for panel b).

Definitions
-----------
Panel a: for each disease, the mean over its genes with score > 0 of the Open Targets
``genetic_association`` data-type score (direct associations). Rare/common labels: MONDO
rare/common table, Orphanet NON RARE -> common and OBSOLETE -> excluded, expert relabel
rounds 1 (common -> rare) and 2 (-> common). Diseases whose name or description contains
'cancer' are excluded. All diseases are plotted (no subsampling) as % of diseases within each group;
the histogram uses fixed bins of width 0.05 on [0, 1].
Panel b: a gene is a rare-disease gene if a rare disease has score > 0.5 for it from a
rare-disease source (Orphanet, ClinVar/EVA, ClinGen, Gene2Phenotype, Genomics England,
UniProt variants, UniProt literature); it is a common-disease gene if a common disease has a
genome-wide association score > 0.5 for it (GWAS credible sets; ``ot_genetics_portal`` in
releases before 24.06). Labels: MONDO rare/common table only, all labelled diseases.
SI panel: share of direct associations contributed by each germline genetic-association
source, among the diseases of panel a.
"""
import argparse
import glob
import json
import os
import re

import numpy as np
import pandas as pd

SCORE = "genetic_association_score"
AVG = "avg_" + SCORE
MAX_DIS_GENE_ASSOC = 30000
BIN_EDGES = np.round(np.linspace(0.0, 1.0, 21), 2)  # fixed 0.05 bins on [0, 1]
STRONG = 0.5

GERMLINE_SOURCES = [  # (datasourceId, display label)
    ("gwas_credible_sets", "GWAS credible sets"),
    ("ot_genetics_portal", "GWAS (Open Targets Genetics)"),
    ("gene_burden", "Gene burden"),
    ("eva", "ClinVar"),
    ("genomics_england", "Genomics England"),
    ("gene2phenotype", "Gene2Phenotype"),
    ("orphanet", "Orphanet"),
    ("clingen", "ClinGen"),
    ("uniprot_literature", "UniProt literature"),
    ("uniprot_variants", "UniProt variants"),
]
RARE_SOURCES = {"orphanet", "gene2phenotype", "clingen", "genomics_england", "eva",
                "uniprot_variants", "uniprot_literature"}
GWAS_SOURCES = {"gwas_credible_sets", "ot_genetics_portal"}


# --------------------------------------------------------------------------- loaders
def _read_parquet_dir(path, columns=None):
    files = sorted(glob.glob(os.path.join(path, "*.parquet")))
    if not files:
        raise FileNotFoundError(f"no parquet files in {path}")
    return pd.concat([pd.read_parquet(f, columns=columns) for f in files], ignore_index=True)


# datatype of each germline source in releases >= 25.x (from the evidence parquet ``datatypeId``);
# 2024 exports carry the datatype column themselves
DATATYPE_OF_SOURCE = {"gwas_credible_sets": "genetic_association", "gene_burden": "genetic_association",
                      "eva": "genetic_association", "orphanet": "genetic_association",
                      "uniprot_variants": "genetic_association", "clingen": "genetic_literature",
                      "gene2phenotype": "genetic_literature", "genomics_england": "genetic_literature",
                      "uniprot_literature": "genetic_literature"}


def _xref_first(xs, prefix):
    if xs is None:
        return None
    if isinstance(xs, str):
        xs = re.findall(r"[A-Za-z0-9_.-]+:[A-Za-z0-9_.-]+", xs)
    for x in xs:
        if isinstance(x, str) and x.startswith(prefix):
            return x
    return None


def load_release(ot_root, fmt, ga_datatypes=("genetic_association",)):
    """Return (diseases, ga, by_source) in one schema for both input formats.

    diseases: id, name, description, MONDO_id, Orphanet_id
    ga:       diseaseid, targetid, genetic_association_score (direct, per data type; when more than
              one data type is given, the per-pair maximum is used)
    by_source: diseaseid, targetid, datatypeid, datasourceid, score
    """
    if fmt == "csv2024":
        d = pd.read_csv(os.path.join(ot_root, "opentargets-diseases.csv"))
        ga = pd.read_csv(os.path.join(ot_root, "ot_assoc_lookup.csv"),
                         usecols=["diseaseid", "targetid", SCORE])
        src = pd.read_csv(os.path.join(ot_root, "ot_gene_association_by_datasource.csv"),
                          usecols=["diseaseid", "targetid", "datatypeid", "datasourceid", "score"])
        return d, ga, src
    d = _read_parquet_dir(os.path.join(ot_root, "disease"), ["id", "name", "description", "dbXRefs"])
    # 2024 rule: an id that is itself MONDO/Orphanet maps to itself, otherwise the first cross-reference
    d["MONDO_id"] = np.where(d["id"].str.startswith("MONDO_"), d["id"],
                             d["dbXRefs"].apply(lambda x: _xref_first(x, "MONDO:")))
    d["Orphanet_id"] = np.where(d["id"].str.startswith("Orphanet_"), d["id"],
                                d["dbXRefs"].apply(lambda x: _xref_first(x, "Orphanet:")))
    d = d.drop(columns=["dbXRefs"])
    a = _read_parquet_dir(os.path.join(ot_root, "association_by_datatype_direct"),
                          ["diseaseId", "targetId", "aggregationType", "aggregationValue", "associationScore"])
    a = a[(a["aggregationType"] == "datatypeId") & a["aggregationValue"].isin(list(ga_datatypes))]
    ga = (a.groupby(["diseaseId", "targetId"], as_index=False)["associationScore"].max()
          .rename(columns={"diseaseId": "diseaseid", "targetId": "targetid", "associationScore": SCORE}))
    s = _read_parquet_dir(os.path.join(ot_root, "association_by_datasource_direct"),
                          ["diseaseId", "targetId", "aggregationType", "aggregationValue", "associationScore"])
    s = s[(s["aggregationType"] == "datasourceId") & s["aggregationValue"].isin(list(DATATYPE_OF_SOURCE))]
    src = pd.DataFrame({"diseaseid": s["diseaseId"].values, "targetid": s["targetId"].values,
                        "datasourceid": s["aggregationValue"].values, "score": s["associationScore"].values})
    src["datatypeid"] = src["datasourceid"].map(DATATYPE_OF_SOURCE)
    return d, ga, src


def check_2024_id_rule(d):
    """Re-derive MONDO_id / Orphanet_id of the 2024 export from dbxrefs with the rule used for later releases."""
    m = np.where(d["id"].str.startswith("MONDO_"), d["id"], d["dbxrefs"].apply(lambda x: _xref_first(x, "MONDO:")))
    o = np.where(d["id"].str.startswith("Orphanet_"), d["id"],
                 d["dbxrefs"].apply(lambda x: _xref_first(x, "Orphanet:")))
    nm = lambda v: pd.Series(v, index=d.index).astype(str).str.replace("_", ":")
    return dict(mondo_match=int((nm(m) == nm(d["MONDO_id"])).sum()), mondo_nonnull=int(d["MONDO_id"].notna().sum()),
                mondo_derived_nonnull=int(pd.Series(m).notna().sum()),
                orphanet_match=int((nm(o) == nm(d["Orphanet_id"])).sum()),
                orphanet_nonnull=int(d["Orphanet_id"].notna().sum()),
                orphanet_derived_nonnull=int(pd.Series(o).notna().sum()))


def _crosswalk(labels_root):
    p = os.path.join(labels_root, "id_crosswalk.csv")
    if not os.path.exists(p):
        return {}
    cw = pd.read_csv(p)
    return dict(zip(cw.iloc[:, 0], cw.iloc[:, 1]))


def _with_successors(ids, cw):
    ids = set(ids)
    return ids | {cw[i] for i in ids if i in cw}


# --------------------------------------------------------------------------- panel a
def build_panel_a_population(d, ga, labels_root):
    """Disease-level table for panel a (2024 logic, unchanged)."""
    ot = d.copy()
    m = ot["MONDO_id"].fillna("").str.contains("_")
    ot.loc[m, "MONDO_id"] = ot.loc[m, "MONDO_id"].str.replace("_", ":")
    m = ot["Orphanet_id"].notna()
    ot.loc[m, "orphanet_disease_id"] = (
        ot.loc[m, "Orphanet_id"].str.replace(":", "_").apply(lambda x: x.split("_")[-1]).astype(int))

    dg = ga[ga["diseaseid"].isin(ot["id"])]
    grp = (dg[dg[SCORE] > 0].groupby("diseaseid")
           .agg(dis_gene_count=("targetid", np.size), **{AVG: (SCORE, "mean")}).reset_index())

    mondo = pd.read_parquet(os.path.join(labels_root, "mondo_rare_common.parquet"))
    mondo["rare"] = mondo["rare"].str.replace("NA", "NA: Cant Assign Using MONDO")
    mondo["MONDO_id"] = mondo["id"].apply(lambda x: x.split("/")[-1].replace("_", ":"))
    dm = pd.merge(ot, mondo.drop(columns=["id"]), on="MONDO_id", how="left")
    dm["rare"] = dm["rare"].fillna("NA: OT missing MONDO")

    idout = pd.read_csv(os.path.join(labels_root, "orphanet_idout.csv"))
    nm = idout["orphanet_disease_name"]
    no_longer_rare = idout.loc[nm.str.contains("NON RARE", case=False, na=False), "orphanet_disease_id"]
    obsolete_rare = idout.loc[nm.str.contains("OBSOLETE", case=False, na=False), "orphanet_disease_id"]
    dm.loc[dm["orphanet_disease_id"].isin(obsolete_rare), "rare"] = "Orpha obsolete"
    dm.loc[dm["orphanet_disease_id"].isin(no_longer_rare), "rare"] = "common"

    g = pd.merge(dm.rename(columns={"id": "diseaseid"}), grp, on="diseaseid", how="left")
    n0 = len(g)
    cancer = g["name"].str.contains("cancer", na=False) | g["description"].astype(str).str.contains("cancer")
    n_cancer = int(cancer.sum())
    g = g[~cancer].copy()
    g[["dis_gene_count", AVG]] = g[["dis_gene_count", AVG]].fillna(0)

    cw = _crosswalk(labels_root)
    r1 = pd.read_excel(os.path.join(labels_root, "relabel_round1.xlsx"))
    r1 = r1[(r1["rare"] == "common") & (r1["rare_reclass"] == "rare")]
    r1_ids = _with_successors(r1["diseaseid"], cw)
    m1 = g["diseaseid"].isin(r1_ids)
    g.loc[m1, "rare"] = "rare"
    r2 = pd.read_excel(os.path.join(labels_root, "relabel_round2.xlsx"), sheet_name="All Data")
    r2_ids = _with_successors(r2.loc[r2["Babak relabel"] == "common", "OT diseaseid"].unique(), cw)
    m2 = g["diseaseid"].isin(r2_ids)
    g.loc[m2, "rare"] = "common"

    keep = (g["dis_gene_count"] < MAX_DIS_GENE_ASSOC) & (g["obsolete"] == 0) & g["rare"].isin(["common", "rare"])
    labels = dict(zip(g.loc[keep, "diseaseid"], g.loc[keep, "rare"]))  # panel b uses the same labels
    out = g[keep & (g[AVG] > 0)][["diseaseid", "name", "rare", "dis_gene_count", AVG]].reset_index(drop=True)
    summ = dict(n_ot_diseases=int(n0), n_cancer_excluded=n_cancer, n_relabel_round1_applied=int(m1.sum()),
                n_relabel_round2_applied=int(m2.sum()),
                n_common=int((out.rare == "common").sum()), n_rare=int((out.rare == "rare").sum()))
    return out, summ, labels


def bin_all(pop):
    """Histogram of all diseases of panel a, per group, as counts and as % of diseases within the group."""
    com = pop[pop.rare == "common"]
    rar = pop[pop.rare == "rare"]
    hc, _ = np.histogram(com[AVG], BIN_EDGES)
    ha, _ = np.histogram(rar[AVG], BIN_EDGES)
    return pd.DataFrame(dict(bin_left=BIN_EDGES[:-1], bin_right=BIN_EDGES[1:], common=hc, rare_all=ha,
                             common_pct=(100 * hc / hc.sum()).round(4), rare_all_pct=(100 * ha / ha.sum()).round(4)))


# --------------------------------------------------------------------------- SI per-source
def source_shares(pop, src):
    rows = []
    lab = dict(GERMLINE_SOURCES)
    s0 = src[src["datasourceid"].isin(lab)]
    for t in ["common", "rare"]:
        s = s0[s0["diseaseid"].isin(pop.loc[pop.rare == t, "diseaseid"])]
        a = (s.groupby(["datasourceid"]).agg(n_associations=("targetid", "size"),
                                              mean_score=("score", "mean"),
                                              datatypes=("datatypeid", lambda x: ";".join(sorted(set(x)))))
             .reset_index())
        a["percent_of_germline_associations"] = (a["n_associations"] / a["n_associations"].sum() * 100).round(2)
        a["label"] = a["datasourceid"].map(lab)
        a["type"] = t
        rows.append(a)
    order = {k: i for i, (k, _) in enumerate(GERMLINE_SOURCES)}
    out = pd.concat(rows)
    out["_o"] = out["datasourceid"].map(order)
    return out.sort_values(["type", "_o"]).drop(columns="_o")


# --------------------------------------------------------------------------- panel b
EXAMPLE_GENES = {"PCSK9": "ENSG00000169174", "LDLR": "ENSG00000130164", "APOB": "ENSG00000084674",
                 "HNF1A": "ENSG00000135100", "GCK": "ENSG00000106633", "HNF4A": "ENSG00000101076",
                 "MC4R": "ENSG00000166603", "LEP": "ENSG00000174697", "POMC": "ENSG00000115138",
                 "TCF7L2": "ENSG00000148737", "LPA": "ENSG00000198670", "ABCG8": "ENSG00000143921",
                 "PPARG": "ENSG00000132170", "KCNJ11": "ENSG00000187486", "ABCC8": "ENSG00000006071",
                 "SLC30A8": "ENSG00000164756", "APOE": "ENSG00000130203", "LPL": "ENSG00000175445",
                 "CFTR": "ENSG00000001626", "HFE": "ENSG00000010704"}


def venn2024_labels(d, mondo_parquet):
    """Labels as in the 2024 shared-gene computation: MONDO table only, no relabels, no cancer filter."""
    mo = pd.read_parquet(mondo_parquet, columns=["id", "obsolete", "rare"])
    mo["m"] = mo["id"].str.split("/").str[-1].str.replace("_", ":")
    mo = mo[mo["obsolete"] == 0]
    lab = dict(zip(mo["m"], mo["rare"]))
    m = np.where(d["id"].str.startswith("MONDO_"), d["id"].str.replace("_", ":"),
                 d["MONDO_id"].astype(object).where(d["MONDO_id"].notna(), None))
    m = pd.Series(m).astype(str).str.replace("_", ":")
    return dict(zip(d["id"], m.map(lab)))


def shared_genes(d, src, L):
    """Genes with score > 0.5 for a rare disease (rare-disease sources) and for a common disease (GWAS)."""
    NM = dict(zip(d["id"], d["name"]))
    s = src[src["datasourceid"].isin(RARE_SOURCES | GWAS_SOURCES)].copy()
    s["grp"] = s["diseaseid"].map(L)
    s = s[s["grp"].isin(["rare", "common"])]
    strong = s[s["score"] > STRONG]
    rg = set(strong.loc[(strong.grp == "rare") & strong.datasourceid.isin(RARE_SOURCES), "targetid"])
    cg = set(strong.loc[(strong.grp == "common") & strong.datasourceid.isin(GWAS_SOURCES), "targetid"])
    counts = dict(n_rare_disease_genes=len(rg), n_common_disease_gwas_genes=len(cg), n_shared=len(rg & cg),
                  n_rare_only=len(rg - cg), n_common_only=len(cg - rg), threshold=STRONG,
                  gwas_sources=sorted(set(strong.datasourceid) & GWAS_SOURCES),
                  datatypes_by_source={k: sorted(v) for k, v in
                                       s.groupby("datasourceid")["datatypeid"].unique().items()})
    ex = strong[strong["targetid"].isin(EXAMPLE_GENES.values())].copy()
    sym = {v: k for k, v in EXAMPLE_GENES.items()}
    ex["gene"] = ex["targetid"].map(sym)
    ex["disease_name"] = ex["diseaseid"].map(NM)
    ex["side"] = np.where((ex.grp == "rare") & ex.datasourceid.isin(RARE_SOURCES), "rare",
                          np.where((ex.grp == "common") & ex.datasourceid.isin(GWAS_SOURCES), "common", ""))
    ex = ex[ex.side != ""].sort_values(["gene", "side", "score"], ascending=[True, True, False])
    genes = pd.DataFrame({"targetid": sorted(rg | cg)})
    genes["rare_disease_gene"] = genes.targetid.isin(rg)
    genes["common_disease_gwas_gene"] = genes.targetid.isin(cg)
    return counts, ex[["gene", "targetid", "side", "diseaseid", "disease_name", "grp", "datasourceid",
                       "datatypeid", "score"]], genes


# --------------------------------------------------------------------------- tables cmd
def cmd_tables(a):
    os.makedirs(a.out, exist_ok=True)
    tag = a.release
    d, ga, src = load_release(a.ot_root, a.input_format, tuple(a.ga_datatypes.split(",")))
    summ0 = dict(ga_datatypes=a.ga_datatypes)
    pop, summ, labels = build_panel_a_population(d, ga, a.labels_root)
    # GWAS-evidence rule: a common disease is kept only if it has at least one direct
    # association with score > 0 from a GWAS-type data source; the rare group is unchanged.
    gw_present = sorted(set(src["datasourceid"].unique()) & GWAS_SOURCES)
    has_gwas = set(src.loc[src["datasourceid"].isin(GWAS_SOURCES) & (src["score"] > 0), "diseaseid"])
    summ0["common_rule"] = "require_gwas" if a.common_require_gwas else "none"
    summ0["gwas_sources_counted"] = gw_present
    if a.common_require_gwas:
        drop = (pop["rare"] == "common") & ~pop["diseaseid"].isin(has_gwas)
        summ0["n_common_before_gwas_rule"] = int((pop["rare"] == "common").sum())
        summ0["n_common_dropped_no_gwas"] = int(drop.sum())
        summ0["n_common_dropped_no_gwas_ge075"] = int((drop & (pop[AVG] >= 0.75)).sum())
        pop.loc[drop, ["diseaseid", "name", AVG, "dis_gene_count"]].to_csv(
            os.path.join(a.out, f"fig1a_common_dropped_no_gwas_OT{tag}.csv"), index=False)
        pop = pop[~drop].reset_index(drop=True)
        labels = {k: v for k, v in labels.items() if v != "common" or k in has_gwas}
        summ["n_common"] = int((pop.rare == "common").sum())
        summ["n_rare"] = int((pop.rare == "rare").sum())
    if a.input_format == "csv2024":
        summ["id_rule_check_2024"] = check_2024_id_rule(d)
    bins = bin_all(pop)
    sub = pop
    sub.to_csv(os.path.join(a.out, f"fig1a_diseases_OT{tag}.csv"), index=False)
    bins.to_csv(os.path.join(a.out, f"fig1a_histogram_bins_OT{tag}.csv"), index=False)
    shares = source_shares(pop, src)
    shares.to_csv(os.path.join(a.out, f"si_source_shares_OT{tag}.csv"), index=False)
    counts, ex, genes = shared_genes(d, src, labels)
    vp = os.path.join(a.labels_root, "mondo_rare_common_venn.parquet")
    if os.path.exists(vp):  # the 2024 shared-gene label definition, for comparison
        counts_v, ex_v, _ = shared_genes(d, src, venn2024_labels(d, vp))
        ex_v.to_csv(os.path.join(a.out, f"fig1b_example_links_venn2024labels_OT{tag}.csv"), index=False)
        summ["panel_b_venn2024_labels"] = counts_v
    ex.to_csv(os.path.join(a.out, f"fig1b_example_links_OT{tag}.csv"), index=False)
    genes.to_csv(os.path.join(a.out, f"fig1b_genes_OT{tag}.csv"), index=False)
    summ.update(summ0)
    summ.update(release=tag, bin_edges=BIN_EDGES.tolist(),
                panel_b=counts)
    json.dump(summ, open(os.path.join(a.out, f"fig1_summary_OT{tag}.json"), "w"), indent=1, default=str)
    print(json.dumps(summ, indent=1, default=str))


# --------------------------------------------------------------------------- figure cmd
BLUE, ORANGE = "#1f77b4", "#ff7f0e"
BLUE_F, ORANGE_F = "#9ecae1", "#fdbe85"


def _style(plt):
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans", "Arial"],
                         "font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8,
                         "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 7,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
                         "pdf.fonttype": 42, "svg.fonttype": "none"})


def _letter(fig, ax, s, dx=-0.06, dy=0.02):
    bb = ax.get_position()
    fig.text(bb.x0 + dx, bb.y1 + dy, s, fontsize=10, fontweight="bold", va="bottom", ha="left")


def _venn_geometry(n1, n2, n12):
    from scipy.optimize import brentq
    r1, r2 = np.sqrt(n1 / np.pi), np.sqrt(n2 / np.pi)

    def lens(dd):
        if dd >= r1 + r2:
            return 0.0
        if dd <= abs(r1 - r2):
            return np.pi * min(r1, r2) ** 2
        a1 = r1 ** 2 * np.arccos((dd ** 2 + r1 ** 2 - r2 ** 2) / (2 * dd * r1))
        a2 = r2 ** 2 * np.arccos((dd ** 2 + r2 ** 2 - r1 ** 2) / (2 * dd * r2))
        a3 = 0.5 * np.sqrt((-dd + r1 + r2) * (dd + r1 - r2) * (dd - r1 + r2) * (dd + r1 + r2))
        return a1 + a2 - a3

    dist = brentq(lambda x: lens(x) - n12, abs(r1 - r2) + 1e-9, r1 + r2 - 1e-9)
    return r1, r2, dist


def _k(n):
    return f"{n:,}"


def cmd_figure(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Circle, FancyBboxPatch, Ellipse
    _style(plt)
    tag = a.release
    T = a.tables
    bins = pd.read_csv(os.path.join(T, f"fig1a_histogram_bins_OT{tag}.csv"))
    summ = json.load(open(os.path.join(T, f"fig1_summary_OT{tag}.json")))
    pb = summ["panel_b"]
    examples = json.load(open(a.examples))  # [{"rare": ..., "genes": [...], "common": ...}, ...]

    fig = plt.figure(figsize=(7.2, 6.0))
    gs = fig.add_gridspec(2, 2, left=0.08, right=0.98, top=0.93, bottom=0.07, wspace=0.28, hspace=0.42,
                          width_ratios=[1, 1.08])

    # ---- a: histogram
    ax = fig.add_subplot(gs[0, 0])
    w = bins.bin_right - bins.bin_left
    ax.bar(bins.bin_left, bins.common_pct, width=w, align="edge", color=BLUE, alpha=0.55, edgecolor=BLUE, lw=0.5,
           label="common")  # H46: no n in panel a
    ax.bar(bins.bin_left, bins.rare_all_pct, width=w, align="edge", color=ORANGE, alpha=0.5, edgecolor=ORANGE,
           lw=0.5, label="rare")  # H46: no n in panel a
    ax.set_xlim(0, 1)
    ax.set_xticks(np.arange(0, 1.01, 0.2))
    ax.set_xlabel("Mean Open Targets genetic association score per disease")
    ax.set_ylabel("Diseases in group (%)")
    ax.legend(frameon=False, loc="upper center", handlelength=1.2)
    ax.set_title("Rare diseases tend to have stronger genetic evidence", loc="left")
    ax.margins(y=0.08)
    _letter(fig, ax, "a")

    # ---- b: Venn + examples
    gsb = gs[0, 1].subgridspec(1, 2, width_ratios=[0.9, 1.1], wspace=0.08)
    axv = fig.add_subplot(gsb[0, 0])
    axe = fig.add_subplot(gsb[0, 1])
    n1, n2, n12 = pb["n_rare_disease_genes"], pb["n_common_disease_gwas_genes"], pb["n_shared"]
    r1, r2, dist = _venn_geometry(n1, n2, n12)
    s = 1.0 / (r1 + r2 + dist)
    c1, c2 = (0.0, 0.0), (dist * s, 0.0)
    axv.add_patch(Circle(c1, r1 * s, fc=ORANGE_F, ec=ORANGE, lw=0.8, alpha=0.75))
    axv.add_patch(Circle(c2, r2 * s, fc=BLUE_F, ec=BLUE, lw=0.8, alpha=0.6))
    xo = (c2[0] - r2 * s + c1[0] + r1 * s) / 2  # centre of the lens
    axv.text((c1[0] - r1 * s + c2[0] - r2 * s) / 2, 0, _k(n1 - n12), ha="center", va="center", fontsize=7)
    axv.text((c1[0] + r1 * s + c2[0] + r2 * s) / 2, 0, _k(n2 - n12), ha="center", va="center", fontsize=7)
    axv.text(xo, 0, _k(n12), ha="center", va="center", fontsize=7, fontweight="bold")
    axv.annotate("shared", (xo, r1 * s * 0.22), (xo, r1 * s + 0.06), ha="center", va="bottom", fontsize=7,
                 fontweight="bold", arrowprops=dict(arrowstyle="-", lw=0.5, color="0.3"))
    xm = (c1[0] - r1 * s + c2[0] + r2 * s) / 2
    axv.text(xm, -r1 * s - 0.05, f"rare-disease genes (n = {_k(n1)})", ha="center", va="top", color=ORANGE, fontsize=6.5)
    axv.text(xm, -r1 * s - 0.15, f"common-disease GWAS genes\n(n = {_k(n2)})", ha="center", va="top", color=BLUE,
             fontsize=6.5)
    axv.set_xlim(c1[0] - r1 * s - 0.05, c2[0] + r2 * s + 0.05)
    axv.set_ylim(-r1 * s - 0.40, r1 * s + 0.16)
    axv.set_aspect("equal")
    axv.axis("off")
    bbc = gs[0, 1].get_position(fig)
    bba = ax.get_position()
    fig.text(bbc.x0, bba.y1 + 0.012, "Rare and common diseases share genes", fontsize=8, va="bottom", ha="left")
    fig.text(bbc.x0 - 0.035, bba.y1 + 0.02, "b", fontsize=10, fontweight="bold", va="bottom", ha="left")

    axe.set_xlim(0, 1)
    axe.set_ylim(0, 1)
    axe.axis("off")
    nex = len(examples)
    hbox, gap = 0.085, 0.075
    block = 2 * hbox + 0.075
    top = 0.5 + (nex * block + (nex - 1) * gap) / 2
    for i, e in enumerate(examples):
        y0 = top - i * (block + gap)
        for (yy, txt, fc, ec) in [(y0 - hbox, e["rare"], ORANGE_F, ORANGE), (y0 - block, e["common"], BLUE_F, BLUE)]:
            axe.add_patch(FancyBboxPatch((0.02, yy), 0.96, hbox, boxstyle="round,pad=0.005,rounding_size=0.02",
                                         fc=fc, ec=ec, lw=0.6, alpha=0.6))
            axe.text(0.5, yy + hbox / 2, txt, ha="center", va="center", fontsize=6)
        axe.text(0.5, y0 - hbox - 0.0375, ", ".join(e["genes"]), ha="center", va="center", fontsize=6.5,
                 style="italic")

    # ---- c: long tail
    axc = fig.add_subplot(gs[1, 0])
    x = np.linspace(0.02, 1, 400)
    y = 1 / (x * 12 + 0.25) ** 1.6
    cut = 0.25
    axc.plot(x, y, color="0.25", lw=1.0)
    axc.fill_between(x, 0, y, where=x <= cut, color=ORANGE_F, alpha=0.8, lw=0)
    axc.fill_between(x, 0, y, where=x >= cut, color=BLUE_F, alpha=0.8, lw=0)
    axc.set_xticks([])
    axc.set_yticks([])
    axc.set_xlabel("Available data per disease \u2192")
    axc.set_ylabel("Number of diseases")
    axc.set_xlim(0, 1)
    axc.set_ylim(0, y.max() * 1.05)
    axc.text(0.12, y.max() * 0.62, "Rare: too little data per disease\nto learn each one alone",
             color="#a6510b", fontsize=7, va="bottom")
    axc.text(0.36, y.max() * 0.2, "Common: data-rich, but models can learn\nshortcuts (popularity, ontology distance)",
             color="#1a5a8a", fontsize=7, va="bottom")
    axc.set_title("Most diseases are data-poor", loc="left")
    _letter(fig, axc, "c")

    # ---- d: gap schematic
    axd = fig.add_subplot(gs[1, 1])
    axd.set_xlim(0, 1)
    axd.set_ylim(0, 1)
    axd.axis("off")
    axd.set_title("What is missing", loc="left")
    axd.text(0.22, 0.92, "Broad similarity\n(genes, phenotypes)", ha="center", va="top", fontweight="bold",
             color="0.3")
    axd.text(0.76, 0.92, "Specific shared\nmechanism (pathway\nor target)", ha="center", va="top",
             fontweight="bold", color="#2e7d4f")
    axd.add_patch(Ellipse((0.25, 0.57), 0.40, 0.34, fc="none", ec="0.35", ls="--", lw=0.8))
    for (xx, yy, col) in [(0.16, 0.62, ORANGE_F), (0.34, 0.64, BLUE_F), (0.25, 0.49, BLUE_F)]:
        axd.add_patch(Ellipse((xx, yy), 0.09, 0.075, fc=col, ec="none"))
    axd.plot([0.6, 0.9], [0.57, 0.57], color="#2e7d4f", lw=1.0, zorder=1)
    axd.add_patch(Ellipse((0.6, 0.57), 0.09, 0.075, fc=ORANGE_F, ec="none", zorder=2))
    axd.add_patch(Ellipse((0.9, 0.57), 0.09, 0.075, fc=BLUE_F, ec="none", zorder=2))
    axd.add_patch(FancyBboxPatch((0.68, 0.535), 0.14, 0.07, boxstyle="round,pad=0.005", fc="#cfe6d8",
                                 ec="#2e7d4f", lw=0.6, zorder=3))
    axd.text(0.75, 0.57, "module", ha="center", va="center", fontsize=7, zorder=4)
    axd.text(0.5, 0.17, "Needed: one representation of rare and common\n"
                        "diseases that (i) captures both broad and specific\n"
                        "shared mechanisms, (ii) is stable across training\n"
                        "runs and (iii) is tested against shortcuts\n"
                        "(popularity, ontology distance).",
             ha="center", va="center", fontsize=6.5,
             bbox=dict(boxstyle="square,pad=0.5", fc="0.96", ec="0.7", lw=0.6))
    _letter(fig, axd, "d", dx=-0.035)

    os.makedirs(a.out, exist_ok=True)
    fig.savefig(os.path.join(a.out, "knowledge_transfer.png"), dpi=300)
    fig.savefig(os.path.join(a.out, "knowledge_transfer.pdf"))
    _check_overlaps(fig, os.path.join(a.out, "knowledge_transfer.overlaps.txt"))
    plt.close(fig)

    # ---- SI per-source panel
    sh = pd.read_csv(os.path.join(T, f"si_source_shares_OT{tag}.csv"))
    labs = [l for _, l in GERMLINE_SOURCES if l in set(sh.label)]
    fig2, ax2 = plt.subplots(figsize=(3.8, 3.0))
    yy = np.arange(len(labs))
    for off, t, col in [(-0.2, "common", BLUE), (0.2, "rare", ORANGE)]:
        v = sh[sh.type == t].set_index("label").reindex(labs)["percent_of_germline_associations"].fillna(0)
        ax2.barh(yy + off, v.values, height=0.4, color=col, label=t)
    ax2.set_yticks(yy)
    ax2.set_yticklabels(labs)
    ax2.invert_yaxis()
    ax2.set_xlabel("Genetic associations (%)")
    ax2.legend(frameon=False, loc="lower right")
    ax2.set_title("Sources of genetic evidence differ\nbetween rare and common diseases", loc="left")
    ax2.set_xlim(0, 100)
    fig2.tight_layout()
    fig2.savefig(os.path.join(a.out, "si_prevalence_sources.png"), dpi=300)
    fig2.savefig(os.path.join(a.out, "si_prevalence_sources.pdf"))
    _check_overlaps(fig2, os.path.join(a.out, "si_prevalence_sources.overlaps.txt"))
    plt.close(fig2)

    # ---- Supplementary Figure S1 composite: a = existing prevalence plot (unchanged image), b = source shares
    if getattr(a, "s1a", None):
        img = plt.imread(a.s1a)
        h, w = img.shape[:2]
        fig3 = plt.figure(figsize=(7.2, 2.75))
        gs3 = fig3.add_gridspec(1, 2, left=0.03, right=0.98, top=0.86, bottom=0.15, wspace=0.42,
                                width_ratios=[1.45, 1])
        axa = fig3.add_subplot(gs3[0, 0])
        axa.imshow(img, interpolation="lanczos")
        axa.set_aspect("equal")
        axa.axis("off")
        axb = fig3.add_subplot(gs3[0, 1])
        for off, t, col in [(-0.2, "common", BLUE), (0.2, "rare", ORANGE)]:
            v = sh[sh.type == t].set_index("label").reindex(labs)["percent_of_germline_associations"].fillna(0)
            axb.barh(yy + off, v.values, height=0.4, color=col, edgecolor="0.15", lw=0.4,
                     label=f"{t} (n = {_k(summ['n_' + t])})")
        axb.set_yticks(yy)
        axb.set_yticklabels(labs, fontsize=6.5)
        axb.tick_params(axis="x", labelsize=6.5)
        axb.invert_yaxis()
        axb.set_xlim(0, 100)
        axb.set_xlabel("Genetic associations (%)", fontsize=7)
        axb.legend(frameon=False, loc="center right", fontsize=6.5, handlelength=1.2)
        axb.set_title("Sources of genetic evidence", loc="left")
        pa, pbb = axa.get_position(), axb.get_position()
        fig3.text(pa.x0, 0.97, "a", fontsize=10, fontweight="bold", va="top", ha="left")
        fig3.text(pbb.x0 - 0.19, 0.97, "b", fontsize=10, fontweight="bold", va="top", ha="left")
        fig3.savefig(os.path.join(a.out, "si_S1_composite.png"), dpi=300)
        fig3.savefig(os.path.join(a.out, "si_S1_composite.pdf"))
        _check_overlaps(fig3, os.path.join(a.out, "si_S1_composite.overlaps.txt"))
        plt.close(fig3)


def _check_overlaps(fig, path):
    import matplotlib as mpl
    r = fig.canvas.get_renderer()
    texts = [(t, t.get_window_extent(r)) for t in fig.findobj(mpl.text.Text)
             if t.get_text().strip() and t.get_visible()]
    bad = []
    for i, (ta, ba) in enumerate(texts):
        if not fig.bbox.contains(ba.x0, ba.y0) or not fig.bbox.contains(ba.x1, ba.y1):
            bad.append(f"outside figure: {ta.get_text()!r}")
        for tb, bb in texts[i + 1:]:
            if ba.overlaps(bb):
                bad.append(f"overlap: {ta.get_text()!r} <> {tb.get_text()!r}")
    with open(path, "w") as f:
        f.write("\n".join(bad) + ("\n" if bad else ""))
    print(f"{os.path.basename(path)}: {len(bad)} issue(s)")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    t = sp.add_parser("tables")
    t.add_argument("--ot-root", required=True)
    t.add_argument("--labels-root", required=True)
    t.add_argument("--release", required=True, help="release tag written into file names, e.g. 26.09")
    t.add_argument("--input-format", choices=["parquet", "csv2024"], default="parquet")
    t.add_argument("--ga-datatypes", default="genetic_association",
                   help="comma-separated data types whose per-pair maximum is the panel a score (parquet input)")
    t.add_argument("--common-require-gwas", action="store_true",
                   help="keep a common disease only if it has a direct GWAS-type association (GWAS_SOURCES)")
    t.add_argument("--out", required=True)
    f = sp.add_parser("figure")
    f.add_argument("--tables", required=True)
    f.add_argument("--release", required=True)
    f.add_argument("--examples", required=True, help="JSON list of {rare, genes, common} example boxes")
    f.add_argument("--s1a", default=None, help="existing Supplementary Figure S1 image used unchanged as panel a")
    f.add_argument("--out", required=True)
    a = p.parse_args()
    cmd_tables(a) if a.cmd == "tables" else cmd_figure(a)
