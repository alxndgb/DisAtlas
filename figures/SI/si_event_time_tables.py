"""SI tables for the event-time holdout, sign rule: positive = DisAtlas better in every column.
dE1 in the source is DisAtlas minus comparator on a lower-is-better percentile, so it is negated here
(advantage = comparator E1 - DisAtlas E1). dE2 is already DisAtlas minus comparator on a higher-is-better share.
Usage: python figures/SI/si_event_time_tables.py"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, re, json, pandas as pd
DT, OUT = DERIVED / "si_event_time" / "deltas.tsv", HERE
D = pd.read_csv(DT, sep="\t")
D["a1"], D["a1lo"], D["a1hi"] = -D.dE1, -D.dE1_hi, -D.dE1_lo
D["a2"], D["a2lo"], D["a2hi"] = D.dE2, D.dE2_lo, D.dE2_hi
def fm(x, nd):
    s = f"{x:+.{nd}f}"; s = "0." + "0" * nd if float(s) == 0 else s
    return s.replace("-", "−")
def cell(r, k, nd): return f"{fm(r[k], nd)} [{fm(r[k + 'lo'], nd)}, {fm(r[k + 'hi'], nd)}]"
SETS = [("pct_p1", "Popularity × weighted-degree groups (co-primary)"), ("pct_s3", "Popularity groups (co-primary)"),
        ("pct_s1", "Weighted-degree groups only"), ("pct_s2", "Matched samples of 20 diseases"),
        ("pct_s4", "All candidate diseases, no grouping"), ("pct_s5", "Combined with popularity")]
SHORT = {"pct_p1": "Popularity × weighted degree", "pct_s3": "Popularity", "pct_s4": "No grouping"}
COMP = [("adamic_adar_w", "Weighted Adamic–Adar"), ("adamic_adar", "Adamic–Adar"), ("popularity", "Disease popularity"),
        ("degree", "MechDist degree"), ("strength", "MechDist weighted degree"), ("gene_jaccard", "Shared genes"),
        ("gobp_simgic", "GO biological-process similarity"), ("hpo_simgic", "HPO phenotype similarity")]
g = lambda sc, an, cp, pool: D[(D.scope == sc) & (D.analysis == an) & (D.comparator == cp) & (D.pool == pool)].iloc[0]
L = ["| Comparison set | Comparator | Median-percentile advantage, test pool | Median-percentile advantage, all diseases | Top-10% advantage, test pool | Top-10% advantage, all diseases |",
     "|:--|:--|--:|--:|--:|--:|"]
T1 = []
for an, al in SETS:
    for i, (cp, cl) in enumerate(COMP):
        t, a = g("pooled", an, cp, "test"), g("pooled", an, cp, "all")
        row = [cell(t, "a1", 2), cell(a, "a1", 2), cell(t, "a2", 3), cell(a, "a2", 3)]
        T1.append((an, cp, row)); L.append(f"| {al if i == 0 else ''} | {cl} | " + " | ".join(row) + " |")
cap1 = (": Event-time holdout, all decision years pooled: advantage of DisAtlas over each comparator in every comparison set and both candidate pools.\n"
        "Median-percentile advantage is the comparator's median percentile minus that of DisAtlas, in percentile points; top-10% advantage is DisAtlas minus the comparator in the share of drugs whose best-ranked new trial falls in the top 10%.\n"
        "Positive values mean DisAtlas ranks new trials better; brackets give 95% bootstrap intervals over drugs.\n"
        "In the combined-with-popularity set each comparator is also combined with popularity, and the popularity row compares with popularity alone. {#tbl:SI_event_time_deltas}")
SUB = [("pooled", "All decision years (main analysis)"), ("pooled_excl_2020_2021", "Excluding 2020 and 2021"),
       ("pooled_single_event", "Drug-years with one new disease"), ("pooled_aa_pos", "Diseases with a shared MechDist neighbour"),
       ("pooled_aa_zero", "Diseases without a shared MechDist neighbour"), ("pooled_covered_hpo_simgic", "Diseases with HPO phenotype annotation"),
       ("pooled_covered_gobp_simgic", "Diseases with GO biological-process annotation"), ("pooled_covered_gene_jaccard", "Diseases with gene annotation")]
L2 = ["| Subset of new trials | Comparator | Comparison set | Median-percentile advantage, test pool | Median-percentile advantage, all diseases |", "|:--|:--|:--|--:|--:|"]
T2 = []
def sub_rows(sc, sl, cp, cl):
    for j, an in enumerate(["pct_p1", "pct_s3", "pct_s4"]):
        row = [cell(g(sc, an, cp, "test"), "a1", 2), cell(g(sc, an, cp, "all"), "a1", 2)]
        T2.append((sc, cp, an, row)); L2.append(f"| {sl if j == 0 else ''} | {cl if j == 0 else ''} | {SHORT[an]} | " + " | ".join(row) + " |")
for sc, sl in SUB: sub_rows(sc, sl, "adamic_adar_w", "Weighted Adamic–Adar")
L2.append("| *Other comparators ahead of DisAtlas in a subset under a co-primary set* | | | | |")
other = []
for sc, sl in SUB[1:]:
    for cp, cl in COMP[1:]:
        q = D[(D.scope == sc) & (D.comparator == cp) & D.analysis.isin(["pct_p1", "pct_s3"])]
        if (q.a1hi < 0).any(): other.append((sc, sl, cp, cl))
for sc, sl, cp, cl in other: sub_rows(sc, sl, cp, cl)
cap2 = (": Event-time holdout: Median-percentile advantage of DisAtlas (comparator minus DisAtlas, percentile points) on subsets of new trials, for the two co-primary comparison sets and for ranking with no grouping.\n"
        "The upper block compares with weighted Adamic–Adar; the lower block lists every other subset and comparator for which the comparator ranks new trials above DisAtlas under either co-primary set in either pool.\n"
        "Positive values mean DisAtlas ranks new trials better; brackets give 95% bootstrap intervals over drugs. {#tbl:SI_event_time_subsets}")
open(f"{OUT}/SI_event_time_tables.md", "w").write("\n".join(L) + "\n\n" + cap1 + "\n\n" + "\n".join(L2) + "\n\n" + cap2 + "\n")
# numbers used in the section prose
Nn = {}
for an in ["pct_p1", "pct_s3", "pct_s4"]:
    for pool in ["test", "all"]:
        q = D[(D.comparator == "adamic_adar_w") & (D.analysis == an) & (D.pool == pool) & D.scope.str.isnumeric()]
        Nn[f"years_{an}_{pool}"] = dict(n=len(q), favour=int((q.a1lo > 0).sum()), against=int((q.a1hi < 0).sum()),
                                        against_years=q[q.a1hi < 0].scope.tolist(), includes0=int(((q.a1lo <= 0) & (q.a1hi >= 0)).sum()))
json.dump(dict(cells=len(T1) + len(T2), other_rows=[o[:3:2] for o in other], **Nn),
          open(f"{OUT}/SI_event_time_tables_numbers.json", "w"), indent=1)
print(json.dumps(dict(cells=len(T1) + len(T2), other=[(o[0], o[2]) for o in other], **Nn)))
