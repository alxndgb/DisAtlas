"""SI tables of the event-time holdout added with the odds-ratio design of Figure 7:
  tbl:SI_nmo_ms_by_slice           eculizumab, DisAtlas rank of MS and NMO by decision year (held-out candidates)
  tbl:SI_event_time_or_yearly      conditional-logit odds ratio per 1 SD, DisAtlas and weighted Adamic-Adar, by decision year
  tbl:SI_event_time_auc            per-drug probability by method, comparison set and candidate pool, pooled
  tbl:SI_event_time_auc_yearly     per-drug probability of DisAtlas by decision year and comparison set
  tbl:SI_event_time_early_years    early decision years (1996, 1997, 1999) and 2000
  tbl:SI_event_time_small_sets     matched-set sizes by decision year
Every cell is read from data/derived/fig7_prospective/ and the numbers quoted in the SI prose are asserted.
Usage: python figures/SI/si_event_time_or_tables.py [OUTDIR]  -> SI_event_time_or_tables.md, SI_event_time_or_tables_numbers.json
"""
import json, os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
import disatlas.deterministic  # noqa: E402,F401
D = os.path.join(ROOT, "data", "derived", "fig7_prospective")
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
tsv = lambda f, **k: pd.read_csv(os.path.join(D, f), sep="\t", **k)
f1 = lambda x: f"{x:.1f}"; f2 = lambda x: f"{x:.2f}"; f3 = lambda x: f"{x:.3f}"
cell3 = lambda t: f"{f3(t[0])} ({f3(t[1])}–{f3(t[2])})"
PLOT_START = 2000
CAP = {
    'SI_nmo_ms_by_slice': ': Eculizumab, DisAtlas rank of MS and NMO among diseases held out of benchmark tuning, by decision year.\nEach decision year uses the map built from evidence up to the year before.\nRank 1 is the closest disease to PNH, and top % is the rank as a share of candidates.\nNot in map means that NMO has no node in that snapshot.\nThere is no row for 1998. {#tbl:SI_nmo_ms_by_slice}',
    'SI_event_time_or_yearly': ': Event-time holdout, odds that a disease not yet trialed with a drug is the one it is trialed in next, per 1 SD of score, each adjusted for popularity, by decision year and pooled, among diseases held out of benchmark tuning.\nDisAtlas values are the dark blue series of [Figure @fig:event_time]b; weighted Adamic–Adar is computed on the same yearly MechDist graph without learning.\nConditional logistic model with one stratum per drug and decision year (Efron ties), popularity as a natural spline plus a zero indicator; 95% intervals from robust standard errors clustered by drug.\nDecision years share drugs, so the yearly estimates are not independent. {#tbl:SI_event_time_or_yearly}',
    'SI_event_time_auc': ': Event-time holdout, all decision years pooled, on the per-drug probability scale: probability that each method scores the next trialed disease above a matched disease not yet trialed with the drug (random ranking, 0.5), averaged over drug–decision-year pairs, and the paired difference from DisAtlas.\nBrackets give 95% intervals from a bootstrap over drugs (2,000 replicates, the replicates of the published analysis).\nThe pre-specified test is the median-percentile comparison with weighted Adamic–Adar ([Supplementary Table @tbl:SI_event_time_deltas]); this table re-expresses it and the other reference scores.\nShared-gene, GO biological-process and HPO phenotype similarity use current annotation releases for every decision year, so unlike the other scores they are not limited to earlier evidence. {#tbl:SI_event_time_auc}',
    'SI_event_time_auc_yearly': ': Event-time holdout, DisAtlas, by decision year and pooled: per-drug probability that the next trialed disease is scored above an alternative disease not yet trialed with the drug, in each comparison set, in diseases held out of benchmark tuning, in the pre-registered matched analysis.\nMean over drugs with its 95% bootstrap interval over drugs (2,000 replicates), and the median and interquartile range of the per-drug values; the last row averages over all drug–decision-year pairs.\nDrugs are counted in the popularity- and connectivity-matched set; a drug whose only new trials of the year had no matched alternative is not counted in that set.\nThere is no row for 1998. {#tbl:SI_event_time_auc_yearly}',
    'SI_event_time_early_years': ": Event-time holdout, early decision years, DisAtlas, diseases held out of benchmark tuning. [Figure @fig:event_time]b plots decision years from 2000 (its first year is shown for comparison); all pooled numbers use the pre-specified decision years 1996, 1997 and 1999–2021.\nTop-10% share, drugs whose best-ranked new trial of the year is in the top 10% of candidates; random, the same share expected under random ranking given each drug's number of new trials and the sizes of its matched groups; probability, the per-drug probability that the new trial is scored above a matched alternative (random, 0.5).\nBrackets give 95% intervals from a bootstrap over drugs (2,000 replicates).\nA trial whose matched group has fewer than 10 candidates cannot reach the top 10% (all years in [Supplementary Table @tbl:SI_event_time_small_sets]).\nThere is no row for 1998.\nNo evidence snapshot ending in 1997 was built, because snapshots before 2000 were built only at the cut-off years of the mechanisms under study and no mechanism had a 1997 cut-off, so decision year 1998 is not scored. {#tbl:SI_event_time_early_years}",
    'SI_event_time_small_sets': ': Event-time holdout, size of the matched sets by decision year, diseases held out of benchmark tuning.\nThe matched set of a new trial is the trialed disease and its matched alternatives not yet trialed with the drug.\nA trial with fewer than 10 candidates cannot reach the top 10% of its ranking; the count is the same in both comparison sets.\nTrials with no alternative are left out of all analyses.\nThe last row covers the decision years plotted in [Figure @fig:event_time]b.\nThere is no row for 1998. {#tbl:SI_event_time_small_sets}',
}
NUM = {}

# ---- tbl:SI_nmo_ms_by_slice ----------------------------------------------------------------------------------------
T = tsv("nmo_ms_by_slice.tsv")
assert sorted(T.y.unique()) == list(range(1996, 2010)), sorted(T.y.unique())
assert list(T[T.note.fillna("") == "no slice"].y) == [1998]
T = T[T.target.notna()]
ms, nmo = T[T.target == "MS"].set_index("y"), T[T.target == "NMO"].set_index("y")
assert ms.embedded.all() and ms.in_test.astype(bool).all()
assert int(nmo[nmo.embedded.astype(bool)].t.min()) == 2006            # NMO has no node before the 2006 snapshot
assert not nmo[nmo.t < 2006].embedded.astype(bool).any() and nmo[nmo.t >= 2006].in_test.astype(bool).all()
pre = ms[ms.t <= 2005]                                                 # every scored decision year before NMO has a node
NUM.update(ms_top_pct_min=float(pre.pct_test.min()), ms_top_pct_max=float(pre.pct_test.max()),
           ms_years_before_first_trial=[2009 - int(pre.index.min()), 2009 - int(pre.index.max())])
assert (f1(pre.pct_test.min()), f1(pre.pct_test.max()), NUM["ms_years_before_first_trial"]) == ("5.0", "21.1", [13, 3])
n07 = nmo.loc[2007]; assert (int(n07.rank_test), int(n07.N_test), f1(n07.pct_test)) == (45, 1728, "2.6")
r_ = lambda v: f"{int(round(v)):,}" if float(v).is_integer() else f"{v:,.1f}"
L = ["| Decision year | Evidence up to | MS rank | MS top % | NMO rank | NMO top % |", "|---:|---:|---:|---:|---:|---:|"]
for y in sorted(ms.index):
    m, n = ms.loc[y], nmo.loc[y]
    nm = (f"{r_(n.rank_test)} of {int(n.N_test):,}", f"{n.pct_test:.1f}") if bool(n.embedded) else ("not in map", "")
    L.append(f"| {y} | {int(m.t)} | {r_(m.rank_test)} of {int(m.N_test):,} | {m.pct_test:.1f} | {nm[0]} | {nm[1]} |")
BLOCKS = [(L, "SI_nmo_ms_by_slice")]

# ---- tbl:SI_event_time_or_yearly -------------------------------------------------------------------------------------
OR = tsv("or_clogit.tsv", dtype={"set": str})
def orv(model, s):
    r = OR[(OR.set == s) & (OR.model == model)]; assert len(r) == 1 and bool(r.converged.iloc[0]), (model, s); r = r.iloc[0]
    return float(r.OR), float(r.lo95), float(r.hi95), int(r.n_cases), int(r.n_drugs)
OYR = list(range(PLOT_START, 2022))
L = ["| Decision year | New trials | Drugs | DisAtlas, odds ratio per 1 SD (95% CI) | Weighted Adamic–Adar, odds ratio per 1 SD (95% CI) |", "|---|---|---|---|---|"]
o2 = lambda t: f"{f2(t[0])} ({f2(t[1])}–{f2(t[2])})"
for y in OYR:
    a, w = orv("M2", str(y)), orv("M2_adamic_adar_w", str(y)); assert a[3:] == w[3:], y
    assert a[1] > w[2], y                                              # DisAtlas interval above weighted Adamic-Adar every year
    L.append(f"| {y} | {a[3]:,} | {a[4]:,} | {o2(a)} | {o2(w)} |")
a, w = orv("M2", "pooled"), orv("M2_adamic_adar_w", "pooled")
L.append(f"| {PLOT_START}–2021 pooled | {a[3]:,} | {a[4]:,} | {o2(a)} | {o2(w)} |")
BLOCKS.append((L, "SI_event_time_or_yearly"))

# ---- tbl:SI_event_time_auc and tbl:SI_event_time_auc_yearly ----------------------------------------------------------
AU = tsv("auc_per_drug.tsv", dtype={"scope": str})
AUr = {(r.pool, r.scope, r.analysis, r.method): (r.auc, r.lo, r.hi, int(r.n_drugs)) for r in AU.itertuples()}
DA = {(r.pool, r.analysis, r.comparator): (r.dAUC, r.lo, r.hi) for r in tsv("auc_deltas.tsv").itertuples()}
BX = {(r.pool, r.analysis, int(r.y)): r for r in tsv("auc_check_yearly_box.tsv").itertuples()}
q = lambda p, a, m="disatlas_ens5", sc="pooled_drug_year": AUr[(p, sc, a, m)]
MLAB = [("disatlas_ens5", "DisAtlas"), ("adamic_adar_w", "Weighted Adamic–Adar (same graph, no learning)"), ("popularity", "Disease popularity"),
        ("gene_jaccard", "Shared genes"), ("gobp_simgic", "GO biological-process similarity"), ("hpo_simgic", "HPO phenotype similarity")]
SETL = {"pct_p1": "Popularity- and connectivity-matched", "pct_s3": "Popularity-matched"}
POOLL = {"test": "Held out of benchmark tuning", "all": "All diseases"}
L = ["| Method | Comparison set | Candidate diseases | Probability (95% CI) | DisAtlas minus method (95% CI) |", "|---|---|---|---|---|"]
for m, lab in MLAB:
    first = True
    for a_ in ("pct_p1", "pct_s3"):
        for p in ("test", "all"):
            d = "—" if m == "disatlas_ens5" else cell3(DA[(p + "|drug_year", a_, m)])
            if m != "disatlas_ens5": assert DA[(p + "|drug_year", a_, m)][1] > 0, (p, a_, m)   # DisAtlas ahead in all four settings
            L.append(f"| {lab if first else ''} | {SETL[a_]} | {POOLL[p]} | {cell3(q(p, a_, m))} | {d} |"); first = False
BLOCKS.append((L, "SI_event_time_auc"))
YRS = sorted(int(s) for (p, s, a_, m) in AUr if p == "test" and a_ == "pct_p1" and m == "disatlas_ens5" and s.isdigit())
assert len(YRS) == 25 and 1998 not in YRS, YRS
YSETS = [("pct_p1", "Popularity- and connectivity-matched"), ("pct_s3", "Popularity-matched"), ("pct_s4", "No matching")]
mq = lambda b: f"{f3(b.median)} [{f3(b.q1)}–{f3(b.q3)}]"
L = ["| Decision year | Drugs | " + " | ".join(f"{l}: mean (95% CI) | {l}: median (interquartile range)" for _, l in YSETS) + " |", "|---|---|" + "---|---|" * len(YSETS)]
for y in YRS:
    cells = []
    for a_, _ in YSETS:
        t = AUr[("test", str(y), a_, "disatlas_ens5")]; b = BX[("test", a_, y)]
        assert abs(b.mean - t[0]) < 1e-9 and int(b.n_drugs) == t[3], (a_, y); assert t[1] > 0.5, (a_, y)
        cells += [cell3(t), mq(b)]
    L.append(f"| {y} | {AUr[('test', str(y), 'pct_p1', 'disatlas_ens5')][3]:,} | " + " | ".join(cells) + " |")
L.append("| All years | " + f"{q('test', 'pct_p1')[3]:,} pairs | " + " | ".join(f"{cell3(q('test', a_))} | —" for a_, _ in YSETS) + " |")
BLOCKS.append((L, "SI_event_time_auc_yearly"))
NUM.update(auc_pooled_p1=q("test", "pct_p1")[:3], auc_pooled_s3=q("test", "pct_s3")[:3], auc_pairs=q("test", "pct_p1")[3])

# ---- tbl:SI_event_time_early_years ------------------------------------------------------------------------------------
EYt = tsv("early_years_metrics.tsv"); EY = {(r.set, int(r.y)): r for r in EYt[EYt.pool == "test"].itertuples()}
E2 = tsv("e2_per_drug.tsv", dtype={"scope": str})
E2 = E2[(E2.pool == "test") & (E2.method == "disatlas_ens5") & (E2.k == 10) & (E2.kind == "e2") & E2.scope.str.isdigit()]
HY = {(r.analysis, int(r.scope)): r.value for r in E2.itertuples()}
for a_ in ("pct_s3", "pct_p1", "pct_s4"):
    for y in (1996, 1997, 1999, 2000): assert abs(EY[(a_, y)].E2 - HY[(a_, y)]) < 5e-4, (a_, y)   # agrees with the SI hits10 figure
EARLY = [1996, 1997, 1999, PLOT_START]
ESET = [("pct_s3", "Popularity-matched"), ("pct_p1", "Popularity- and connectivity-matched"), ("pct_s4", "No matching")]
p3 = lambda v: f"{100 * v:.1f}"
L = ["| Decision year | Evidence up to | Drugs | New trials | New trials per drug | Comparison set | Median candidate diseases | "
     "Trials with fewer than 10 candidates (%) | Drugs with a new trial in the top 10% (%, 95% CI) | Random (%) | Probability (95% CI) |",
     "|---|---|---|---|---|---|---|---|---|---|---|"]
for y in EARLY:
    for k, (a_, lab) in enumerate(ESET):
        r = EY[(a_, y)]
        h = f"| {y} | {int(r.slice[1:])} | {int(r.n_drugs):,} | {int(r.n_events):,} | {r.events_per_drug:.2f} " if k == 0 else "| | | | | "
        L.append(h + f"| {lab} | {float(r.N_median):g} | {p3(r.share_N_lt10)} | {p3(r.E2)} ({p3(r.E2_lo)}–{p3(r.E2_hi)}) | "
                 f"{p3(r.chance)} | {f3(r.auc)} ({f3(r.auc_lo)}–{f3(r.auc_hi)}) |")
BLOCKS.append((L, "SI_event_time_early_years"))

# ---- tbl:SI_event_time_small_sets ------------------------------------------------------------------------------------
SS = tsv("small_sets_by_year.tsv")
assert list(SS.y) == YRS and (SS.p1_lt10 == SS.s3_lt10).all()          # fewer than 10 candidates: same count in both matched sets
for r in SS[SS.y >= PLOT_START].itertuples():
    assert abs(r.p1_lt10_pct - 100 * EY[("pct_p1", r.y)].share_N_lt10) < 0.06 if ("pct_p1", r.y) in EY else True
S20 = SS[SS.y >= PLOT_START]; NT20, NL20 = int(S20.trials.sum()), int(S20.p1_lt10.sum())
assert (NT20, NL20) == (34592, 297) and (S20.p1_noalt == 0).all() and (S20.s3_noalt == 0).all(), (NT20, NL20)
assert NT20 == orv("M2", "pooled")[3]                                 # same new trials as the pooled odds ratio
L = ["| Decision year | New trials | Median matched-set size, popularity- and connectivity-matched | Median matched-set size, popularity-matched | "
     "Trials with fewer than 10 candidates (n, %) | Trials with no alternative (n) |", "|---|---|---|---|---|---|"]
for r in SS.itertuples():
    L.append(f"| {r.y} | {int(r.trials):,} | {float(r.p1_medN):g} | {float(r.s3_medN):g} | {int(r.p1_lt10):,} ({r.p1_lt10_pct:.1f}) | {r.p1_noalt} |")
L.append(f"| {PLOT_START}–{YRS[-1]} | {NT20:,} | — | — | {NL20:,} ({100 * NL20 / NT20:.1f}) | 0 |")
BLOCKS.append((L, "SI_event_time_small_sets"))

# ---- numbers quoted in the SI prose of this section (each list is checked against its manuscript sentence) ----------
ci2 = lambda t: f"{f2(t[0])} (95% CI {f2(t[1])}–{f2(t[2])})"; ci1 = lambda t: f"{f1(t[0])} (95% CI {f1(t[1])}–{f1(t[2])})"
ci3 = lambda t: f"{f3(t[0])} (95% CI {f3(t[1])}–{f3(t[2])})"
P = {}
yw = {y: orv("M2_adamic_adar_w", str(y)) for y in OYR}; ymin = min(yw, key=lambda y: yw[y][0])
assert ymin == 2020 and min(OYR, key=lambda y: orv("M2", str(y))[0]) == 2020
P["waa_or"] = [ci2(w), f"each of the {len(OYR)} decision years"]
P["waa_2020"] = [f"lowest in {ymin} ({f2(yw[ymin][0])}, 95% CI {f2(yw[ymin][1])}–{f2(yw[ymin][2])})"]
AP = tsv("abs_percentiles.tsv", dtype={"scope": str}); AP = AP[(AP.pool == "test") & (AP.scope == "pooled") & (AP.n_rep > 0)]
ap = lambda a_, m: tuple(AP[(AP.analysis == a_) & (AP.method == m)][["E1", "lo", "hi"]].iloc[0])
P["prereg_abs"] = [ci1(ap("pct_p1", "adamic_adar_w")), f"against {f1(ap('pct_p1', 'disatlas_ens5')[0])} for DisAtlas",
                   ci1(ap("pct_s3", "adamic_adar_w")), f"against {f1(ap('pct_s3', 'disatlas_ens5')[0])} among"]
DEC = json.load(open(os.path.join(D, "decision_event_time.json"))); assert DEC["verdict"] == "SUPPORTED" and DEC["n_pass"] == 4
P["prereg_rule"] = [{4: "all four comparisons met it"}[DEC["n_pass"]]]
DL = pd.read_csv(os.path.join(ROOT, "data", "derived", "si_event_time", "deltas.tsv"), sep="\t", dtype={"scope": str})
dl = lambda sc, p, a_, m: DL[(DL.scope == sc) & (DL.pool == p) & (DL.analysis == a_) & (DL.comparator == m)].iloc[0]
adv = lambda p, a_: (-dl("pooled", p, a_, "adamic_adar_w").dE1, -dl("pooled", p, a_, "adamic_adar_w").dE1_hi, -dl("pooled", p, a_, "adamic_adar_w").dE1_lo)
P["prereg_adv"] = [f"{f1(adv('test', 'pct_p1')[0])} percentile points (95% CI {f1(adv('test', 'pct_p1')[1])}–{f1(adv('test', 'pct_p1')[2])})",
                   ci1(adv("test", "pct_s3")), ci1(adv("all", "pct_p1")), ci1(adv("all", "pct_s3")).replace("95% CI ", "95% CI ")]
YY = DL[(DL.pool == "test") & (DL.analysis == "pct_p1") & (DL.comparator == "adamic_adar_w") & DL.scope.str.isdigit()]
P["prereg_years"] = [f"{int((YY.dE1_hi < 0).sum())} of the {len(YY)} decision years", "below zero in " + {0: "none"}[int((YY.dE1_lo > 0).sum())]]
u4 = dl("pooled", "test", "pct_s4", "popularity"); assert u4.comp_E1 < u4.disatlas_E1
P["unmatched"] = [f"median percentile of {f1(u4.comp_E1)}", f"against {f1(u4.disatlas_E1)} for DisAtlas"]
DR = tsv("auc_dropped.tsv"); DR = DR[(DR.pool == "test") & (DR.method == "disatlas_ens5")].set_index("analysis")
assert DR.loc["pct_p1", "n_events_no_negative"] == DR.loc["pct_s3", "n_events_no_negative"]
P["auc_dropped"] = [f"the {int(DR.loc['pct_s3', 'n_events_no_negative'])} new trials without any matched candidate"]
CK = {(r.pool, r.analysis, r.group): (int(r.n), r.mean_auc) for r in tsv("auc_check.tsv").itertuples()}
assert CK[("test", "pct_p1", "per_drug_year")][0] == q("test", "pct_p1")[3] and abs(CK[("test", "pct_p1", "per_drug_year")][1] - q("test", "pct_p1")[0]) < 1e-9
P["auc_pooled"] = [f"all {q('test', 'pct_p1')[3]:,} drug–decision-year pairs", ci3(q("test", "pct_p1")), ci3(q("test", "pct_s3"))]
pd3 = q("test", "pct_s3", sc="pooled"); one, five = CK[("test", "pct_s3", "drugs 1 year")], CK[("test", "pct_s3", "drugs 5+ years")]
assert one[1] > pd3[0] > five[1]
P["auc_per_drug"] = [f"each of the {pd3[3]:,} drugs", ci3(pd3), f"(mean {f3(one[1])}, {one[0]:,} drugs)", f"(mean {f3(five[1])}, {five[0]:,} drugs)"]
gj = q("test", "pct_p1", "gene_jaccard")
assert all(q("test", "pct_p1", m)[0] < gj[0] for m in ("adamic_adar_w", "gobp_simgic", "hpo_simgic", "popularity"))
P["auc_pop"] = [f"reaches {f3(q('test', 'pct_s3', 'popularity')[0])} (popularity-matched) and {f3(q('test', 'pct_p1', 'popularity')[0])}",
                f"at {f3(gj[0])} against {f3(q('test', 'pct_p1')[0])} for DisAtlas"]
g2t, g2a = dl("pooled", "test", "pct_p1", "gene_jaccard"), dl("pooled", "all", "pct_p1", "gene_jaccard")
assert g2t.dE2_lo < 0 < g2t.dE2_hi and g2a.dE2_hi < 0                  # 'interval includes zero ... slightly ahead among all diseases'
P["e2_gene"] = ["interval includes zero", "slightly ahead among all diseases"]
P["ms_proxy"] = [f"from {int(pre.index.min())} to {int(pre.index.max())}", f"{NUM['ms_years_before_first_trial'][0]} to {NUM['ms_years_before_first_trial'][1]} years",
                 f"top {f1(pre.pct_test.min())}–{f1(pre.pct_test.max())}%"]
ER = tsv("ecu_nmo_trial_records.tsv", dtype={"first_year_rule": str})
assert ER[(ER.kind == "nct_present") & (ER.nct_id == "NCT00904826")].first_year_rule.iloc[0] == "0"   # 2009 trial absent from the trial table
first = ER[ER.kind == "record"].sort_values("first_year_rule").iloc[0]
P["ecu_dating"] = [f"records this pair in {first.first_year_rule}", "(NCT00904826)", f"phase {first.phase.replace('III', '3')} trial {first.nct_id}", f"its {first.first_year_rule} {first.source.capitalize()} validation"]
CJ = json.load(open(os.path.join(D, "nmo_case_clogit.json"))); H = json.load(open(os.path.join(D, "nmo_hops.json")))["stripped_training_graph"]
assert H["direct_edge"] is False
ordn = lambda n: f"{n:,}" + ("th" if 10 <= n % 100 <= 20 else {1: "st", 2: "nd", 3: "rd"}.get(n % 10, "th"))
P["nmo_2007"] = [f"In decision year {CJ['year']}, NMO ranks {ordn(int(CJ['disatlas_rank']))} of {int(CJ['N']):,} candidates", f"{ordn(int(CJ['popularity_rank']))} by popularity alone",
                 f"the {CJ['t']} model", f"({H['n_nodes']:,} diseases, {H['n_edges']:,} links) and are {H['hops']} steps apart"]
CS = tsv("nmo_candidate_sets.tsv"); CS = CS[CS.pool == "held-out pool"].set_index("set")
h2, h3 = CS.loc["<=2 hops of PNH"], CS.loc["<=3 hops of PNH"]; AL = tsv("nmo_candidate_sets.tsv"); AL = AL[AL.pool != "held-out pool"]
AL = AL[AL.set == "any distance (current)"].iloc[0]
assert h2.disatlas_rank == h3.disatlas_rank == CJ["disatlas_rank"]
P["nmo_hops"] = [f"leaves the rank at {int(h2.disatlas_rank)} (top {f1(h2.disatlas_top_pct)}% of {int(h2.N):,} and top {f1(h3.disatlas_top_pct)}% of {int(h3.N):,})",
                 f"among all {int(AL.N):,} embedded diseases NMO ranks {ordn(int(AL.disatlas_rank))} (top {f1(100 * AL.disatlas_rank / AL.N)}%)"]
NUM["prose"] = P

out = "\n\n".join("\n".join(rows) + "\n\n" + CAP[lab] for rows, lab in BLOCKS) + "\n"
open(os.path.join(OUT, "SI_event_time_or_tables.md"), "w", encoding="utf-8").write(out)
json.dump(NUM, open(os.path.join(OUT, "SI_event_time_or_tables_numbers.json"), "w"), indent=1)
print("wrote", len(BLOCKS), "tables;", {k: v for k, v in NUM.items()})
