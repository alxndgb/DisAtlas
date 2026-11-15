"""Numbers quoted in the SI event-time section (design, comparison sets, seeds, subsets, years, census, coverage).

Reads the event-time result tables in data/derived/si_event_time and writes SI_event_time_numbers.json next to this
script. Sign rule as in the SI tables: positive = DisAtlas better. For the median percentile the advantage is the
comparator's value minus DisAtlas's (adv = -dE1, interval [-dE1_hi, -dE1_lo]); for the top-10% share it is DisAtlas minus
the comparator (adv = -dE2 as stored comparator-minus-DisAtlas is converted the same way as in si_event_time_tables.py).
Usage: python figures/SI/si_event_time_numbers.py
"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, write_numbers  # noqa: E402
import pandas as pd  # noqa: E402

T = DERIVED / "si_event_time"
D = pd.read_csv(T / "deltas.tsv", sep="\t", dtype={"scope": str})
P = pd.concat([pd.read_csv(T / f"point_estimates_{p}.tsv", sep="\t", dtype={"scope": str}) for p in ("test", "all")])
CEN = pd.read_csv(T / "merge_and_matching_census.tsv", sep="\t", dtype={"scope": str})
COV = pd.read_csv(T / "coverage.tsv", sep="\t")
PR = pd.read_csv(T / "table_pooled_top10_precision.tsv", sep="\t")
Q1 = pd.read_csv(T / "posthoc_covid_2020_counts.tsv", sep="\t")
Q3 = pd.read_csv(T / "posthoc_covid_2020_dE1.tsv", sep="\t")
POOLS = ("test", "all")
r2 = lambda x: round(float(x), 2)
r3 = lambda x: round(float(x), 3)


def adv(analysis, comparator, pool, scope="pooled", endpoint="E1"):
    q = D[(D.analysis == analysis) & (D.comparator == comparator) & (D.pool == pool) & (D.scope == scope)]
    assert len(q) == 1, (analysis, comparator, pool, scope, len(q)); q = q.iloc[0]
    if endpoint == "E1":
        return [r2(-q.dE1), r2(-q.dE1_hi), r2(-q.dE1_lo)]
    return [r3(q.dE2), r3(q.dE2_lo), r3(q.dE2_hi)]


def pe(analysis, method, pool, scope="pooled", col="E1"):
    q = P[(P.analysis == analysis) & (P.method == method) & (P.pool == pool) & (P.scope == scope)]
    assert len(q) == 1, (analysis, method, pool, scope, len(q)); return float(q.iloc[0][col])


N = {}
for pool in POOLS:
    c = CEN[(CEN.pool == pool) & (CEN.scope == "pooled")].iloc[0]
    N[f"design_{pool}"] = dict(n_events=int(c.n_events), n_drugs=int(c.n_drugs))
    N[f"census_{pool}"] = dict(p1_merge0=int(c.p1_merge0), p1_merge1=int(c.p1_merge1), p1_merge2=int(c.p1_merge2), p1_small=int(c.p1_small),
                               s2_widen=int(c.s2_widen), s2_short=int(c.s2_short), s4_n_cand_median=int(c.n_cand_median))
    N[f"census_{pool}"]["s2_drawsd_median_disatlas"] = r2(pe("pct_s2_drawsd_median", "disatlas_ens5", pool))
    for y in ("1996", "1997", "2010"):
        N[f"n_events_{pool}_{y}"] = int(CEN[(CEN.pool == pool) & (CEN.scope == y)].n_events.iloc[0])
    # comparison sets
    N[f"s4_{pool}"] = dict(vs_popularity=adv("pct_s4", "popularity", pool), vs_wAA=adv("pct_s4", "adamic_adar_w", pool),
                           E1_popularity=r2(pe("pct_s4", "popularity", pool)), E1_disatlas=r2(pe("pct_s4", "disatlas_ens5", pool)),
                           E1_wAA=r2(pe("pct_s4", "adamic_adar_w", pool)))
    N[f"s2_{pool}_vs_wAA"] = adv("pct_s2", "adamic_adar_w", pool)
    N[f"s5_{pool}_vs_popularity"] = dict(E1=adv("pct_s5", "popularity", pool), E2=adv("pct_s5", "popularity", pool, endpoint="E2"))
    N[f"p1_{pool}_top10_vs_wAA"] = adv("pct_p1", "adamic_adar_w", pool, endpoint="E2")
    for an in ("prec_s4", "prec_s3"):
        q = PR[(PR.pool == pool) & (PR.analysis == an)].set_index("comparator")
        N[f"{an}_{pool}_mean_pct"] = dict(disatlas=r2(100 * q.DisAtlas_mean.iloc[0]), popularity=r2(100 * q.loc["popularity", "comp_mean"]),
                                          wAA=r2(100 * q.loc["adamic_adar_w", "comp_mean"]))
    # seeds
    for an in ("pct_p1", "pct_s3", "pct_s4"):
        s = [pe(an, f"disatlas_s{k}", pool) for k in (4707, 4708, 4709, 4710, 4711)]
        N[f"seeds_{an}_{pool}"] = dict(seed_min=r2(min(s)), seed_max=r2(max(s)), seed_mean=r2(sum(s) / 5), ens5=r2(pe(an, "disatlas_ens5", pool)),
                                       wAA=r2(pe(an, "adamic_adar_w", pool)), popularity=r2(pe(an, "popularity", pool)))
    # subsets (upper block of the subsets table plus prose)
    for sc in ("pooled_excl_2020_2021", "pooled_single_event", "pooled_aa_pos", "pooled_aa_zero", "pooled_covered_hpo_simgic"):
        N[f"subset_{sc}_{pool}"] = {an: adv(an, "adamic_adar_w", pool, sc) for an in ("pct_p1", "pct_s3", "pct_s4")}
    N[f"subset_single_event_{pool}_AA_p1"] = adv("pct_p1", "adamic_adar", pool, "pooled_single_event")
    N[f"subset_aa_zero_{pool}_p1"] = {cmp: adv("pct_p1", cmp, pool, "pooled_aa_zero") for cmp in ("popularity", "gene_jaccard", "hpo_simgic")}
    N[f"subset_hpo_{pool}_vs_hpo"] = {an: adv(an, "hpo_simgic", pool, "pooled_covered_hpo_simgic") for an in ("pct_p1", "pct_s3", "pct_s4")}
    # individual years
    for an in ("pct_p1", "pct_s3"):
        q = D[(D.comparator == "adamic_adar_w") & (D.analysis == an) & (D.pool == pool) & D.scope.str.isnumeric()]
        lo, hi = -q.dE1_hi, -q.dE1_lo
        N[f"years_{an}_{pool}"] = dict(n=len(q), favour=int((lo > 0).sum()), against=int((hi < 0).sum()), includes0=int(((lo <= 0) & (hi >= 0)).sum()),
                                       against_years=q[hi < 0].scope.tolist())
        for y in ("2020", "1996", "1997"):
            N[f"year_{y}_{an}_{pool}"] = adv(an, "adamic_adar_w", pool, y)
    # coverage
    for m in ("gene_jaccard", "gobp_simgic", "hpo_simgic"):
        N[f"coverage_{m}_{pool}_pct"] = round(100 * float(COV[(COV.pool == pool) & (COV.method == m)].frac_covered.iloc[0]), 1)
q1 = Q1[Q1.pool == "all"].iloc[0]
N["posthoc_2020_all"] = dict(n_events=int(q1.n_events), n_events_covid=int(q1.n_events_covid))
for an in ("pct_p1", "pct_s3"):
    for sc in ("full", "excl_covid_diseases_narrow"):
        q = Q3[(Q3.pool == "all") & (Q3.scenario == sc) & (Q3.analysis == an)]; assert len(q) == 1, (an, sc); q = q.iloc[0]
        N[f"posthoc_2020_all_{sc}_{an}"] = [r2(-q.dE1), r2(-q.dE1_hi), r2(-q.dE1_lo)]
write_numbers(__file__, N, "SI_event_time")
print(len(N), "number groups")
