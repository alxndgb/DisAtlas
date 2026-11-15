"""Quantities for tests/si_text_numbers.tsv (Supplementary Information number registry).

Keys are
  - every key of the main registry (tests/main_quantities.py + value tables of tests/test_main_parity.py), and
  - table cells of the SI / derived tables:  "<relpath>#<row label>|<column>", where <row label> joins the
    non-numeric cells of the row with "/" (suffix "~n" for the n-th repeat of a label), so keys do not depend on
    row order. JSON files are already covered by the main loader as "<relpath without .json>:<dotted.path>".
"""
import csv, glob, math, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

TABLE_GLOBS = ["figures/SI/*.tsv", "tables/SI/*.tsv", "data/derived/**/*.tsv", "data/derived/**/*.csv",
               "figures/main/*.tsv", "tables/main/*.tsv"]
MAX_ROWS = 20000
DROP_LABEL = {"source", "sources", "bootstrap", "inputs", "input", "note", "notes", "path", "file", "files", "resampling", "src"}


def _num(x):
    try:
        v = float(str(x).replace(",", "").replace("\u2212", "-").rstrip("%"))
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def table_quantities(root=ROOT):
    q = {}
    files = sorted({g for p in TABLE_GLOBS for g in glob.glob(os.path.join(root, p), recursive=True)})
    for f in files:
        if os.path.getsize(f) > 8e6:
            continue
        rel = os.path.relpath(f, root)
        with open(f, newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t" if f.endswith(".tsv") else ","))
        rows = [r for r in rows if r and not r[0].startswith("#")]
        if len(rows) < 2 or len(rows) > MAX_ROWS + 1:
            continue
        hdr = rows[0]; seen = {}
        labcols = [i for i, h in enumerate(hdr) if h.strip().lower() not in DROP_LABEL]
        for r in rows[1:]:
            lab = "/".join(r[i] for i in labcols if i < len(r) and r[i] != "" and _num(r[i]) is None and len(r[i]) <= 60)
            k = seen.get(lab, 0); seen[lab] = k + 1
            lab2 = lab if k == 0 else f"{lab}~{k}"
            for c, x in zip(hdr, r):
                v = _num(x)
                if v is not None:
                    q[f"{rel}#{lab2}|{c}"] = v
    return q


def all_si_quantities():
    import main_quantities as MQ
    import test_main_parity as TMP
    q = MQ.all_quantities(); q.update(TMP._value_tables()); q.update(table_quantities())
    return q


DERIVED_OPS = {"min": min, "max": max, "sum": sum, "diff": lambda v: v[0] - v[1], "ratio": lambda v: v[0] / v[1],
               "range": lambda v: max(v) - min(v), "complement": lambda v: 1 - v[0], "absdiff": lambda v: abs(v[0] - v[1]),
               "mean": lambda v: sum(v) / len(v), "npos": lambda v: float(sum(x > 0 for x in v)),
               "neglog10max": lambda v: float("inf") if max(v) == 0 else -math.log10(max(v)),
               "neg": lambda v: -v[0], "count": lambda v: float(len(v)),
               "twice_ratio": lambda v: 2 * v[0] / v[1],  # pairs per disease: each pair counts at both diseases
               "pct100": lambda v: 100 * v[0],  # fraction printed as a bare percentage (no % sign)
               "millions": lambda v: v[0] / 1e6,
               "sum": lambda v: float(sum(v)),
               "others": lambda v: v[0] - 1}  # all diseases except the one itself


def derived(q, path=os.path.join(HERE, "si_derived_quantities.tsv")):
    """Quantities that are a min / max / sum / difference / ratio of other keys (one row per name: name, op, keys
    joined by ' ;; ')."""
    out = {}; q = dict(q)  # later rows may use earlier derived names
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            ks = [k.strip() for k in r["keys"].split(" ;; ")]
            miss = [k for k in ks if k not in q]
            assert not miss, (r["name"], miss)
            out[r["name"]] = q[r["name"]] = DERIVED_OPS[r["op"]]([q[k] for k in ks])
    return out


def computed_quantities(root=ROOT):
    """Quantities the SI prints with no shipped value table, recomputed here from shipped per-pair inputs (point
    estimates only).  Keys start with "CQ:"."""
    import numpy as np
    import pandas as pd
    from scipy.stats import rankdata
    dv = os.path.join(root, "data", "derived", "ixidn_robustness")
    p = pd.read_csv(os.path.join(dv, "ixidn4364_per_pair_scores.tsv.gz"), sep="\t", usecols=["blind_id", "y", "shared_mondo_organsys"])
    s = pd.read_csv(os.path.join(dv, "ixidn4364_semantic_similarity.tsv"), sep="\t")
    d = p.merge(s, on="blind_id", how="left", validate="1:1")
    d = d[d.shared_mondo_organsys.fillna(0).astype(str).str.lower().isin(["1", "1.0", "true"])]

    def auc(y, x):
        y = np.asarray(y); r = rankdata(x); n1 = y.sum(); n0 = len(y) - n1
        return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    out = {"CQ:s10_within_organ_n_pairs": float(len(d))}
    for c, nm in [("hpo_resnik_bma", "hpo"), ("go_resnik_genebma", "go")]:
        x = d[c].fillna(d[c].min() - 1)          # pairs without the annotation rank last (unscorable-last convention)
        out[f"CQ:s10_within_organ_{nm}_auroc"] = auc(d.y, x)

    # S5 degree-distribution aggregates (per-disease degrees; shipped degree_per_disease.tsv / degree_by_channel.tsv)
    dd = pd.read_csv(os.path.join(root, "data", "derived", "degree_distributions", "degree_per_disease.tsv"), sep="\t")
    com = dd.disease_class == "common"
    out["CQ:s5_pct_common"] = float(com.mean() * 100)
    for g in ["broad", "specific", "combined"]:
        top = dd[f"deg_{g}"].sort_values(ascending=False, kind="stable").index[:130]
        out[f"CQ:s5_top130_common_{g}"] = float(com.loc[top].sum())
    out["CQ:s5_pct_deg300_common"] = float((dd.deg_combined[com] >= 300).mean() * 100)
    out["CQ:s5_pct_deg300_rare"] = float((dd.deg_combined[~com] >= 300).mean() * 100)
    dc = pd.read_csv(os.path.join(root, "data", "derived", "mechanism_channels", "degree_by_channel.tsv"), sep="\t", index_col=0)
    out["CQ:s5_spec_deg50"] = float((dc.spec == 50).sum()); out["CQ:s5_spec_deg1_49"] = float(dc.spec.between(1, 49).sum())
    out["CQ:s5_broad_deg3"] = float((dc.broad == 3).sum()); out["CQ:s5_broad_deg1_2"] = float(dc.broad.isin([1, 2]).sum())

    # S14 odds-ratio analysis: weighted Adamic-Adar yearly conditional-logit fits (shipped or_clogit.tsv, model
    # M2_adamic_adar_w) against DisAtlas (M2), decision years 2000-2021
    oc = pd.read_csv(os.path.join(root, "data", "derived", "fig7_prospective", "or_clogit.tsv"), sep="\t", dtype={"set": str})

    def yearly(m):
        y = oc[(oc.model == m) & oc.set.str.isdigit()].assign(y=lambda d: d.set.astype(int))
        return y[(y.y >= 2000) & (y.y <= 2021)].set_index("y").sort_index()
    w, dis = yearly("M2_adamic_adar_w"), yearly("M2")
    assert list(w.index) == list(dis.index) == list(range(2000, 2022))
    out["CQ:s14_waa_or_n_years"] = float(len(w))
    out["CQ:s14_waa_or_min_year"] = float(w.OR.astype(float).idxmin())
    out["CQ:s14_waa_or_years_below_disatlas"] = float((w.hi95.astype(float) < dis.lo95.astype(float)).sum())
    for y in w.index:
        out[f"CQ:s14_waa_or_{y}"], out[f"CQ:s14_waa_or_{y}_lo"], out[f"CQ:s14_waa_or_{y}_hi"] = (
            float(w.loc[y].OR), float(w.loc[y].lo95), float(w.loc[y].hi95))
    return out


def all_si_quantities_with_derived():
    q = all_si_quantities(); q.update(computed_quantities()); q.update(derived(q)); return q
