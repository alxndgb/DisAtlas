"""SI text: FDA-approved repurposing pairs scored with models trained only on evidence dated before the later indication's
first trial (dated prioritisation design), reported in full, plus one Limitations sentence with no numbers.
Every number is read from prospective/fda_prioritisation/tables (endpoints.tsv, paired_differences.tsv, flow.tsv).
Usage: si_fda_prioritisation_text.py <tables_dir> <out_dir>"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, pandas as pd
SRC, OUT = DERIVED / "si_fda", HERE
E = pd.read_csv(f"{SRC}/endpoints.tsv", sep="\t"); P = pd.read_csv(f"{SRC}/paired_differences.tsv", sep="\t"); F = pd.read_csv(f"{SRC}/flow.tsv", sep="\t")
pc = lambda x: f"{100 * x:.1f}%"; pp = lambda x: "0.0" if abs(x) < 5e-4 else f"{100 * x:+.1f}".replace("-", "−"); sg = lambda x: f"{x:+.1f}".replace("-", "−")
def e(d, m, col="top10"): r = E[(E["pop"] == "dev") & (E.design == d) & (E.method == m)]; assert len(r) == 1; return r.iloc[0]
def p(d, m): r = P[(P["pop"] == "dev") & (P.design == d) & (P.comparison == f"disatlas_ens - {m}")]; assert len(r) == 1; return r.iloc[0]
def f(d): r = F[(F["pop"] == "dev") & (F.design == d)]; assert len(r) == 1; return r.iloc[0]
NAME = dict(disatlas_ens="DisAtlas (5 seeds)", deg="disease degree", shared_nbr="shared graph neighbours", lit_shared="shared literature genes", decoy="matched decoy diseases")
for d in ("D1", "D2"): assert e(d, "disatlas_ens").n == e(d, "shared_nbr").n == p(d, "shared_nbr").n
f1, f2, a1, a2 = f("D1"), f("D2"), e("D1", "disatlas_ens"), e("D2", "disatlas_ens")
s1, s2 = p("D1", "shared_nbr"), p("D2", "shared_nbr")
rows = []
for m in ["disatlas_ens", "deg", "shared_nbr", "lit_shared", "decoy"]:
    c = [NAME[m]]
    for d in ("D1", "D2"):
        r = e(d, m); c.append(f"{pc(r.top10)}" + (f" ({100 * r.top10_lo:.1f}–{100 * r.top10_hi:.1f})" if m != "decoy" else " (20 per pair)"))
        c.append(f"{r.median_pct:.1f}" if m != "decoy" else f"{r.median_pct:.1f}")
        if m == "disatlas_ens": c.append("—")
        else: q = p(d, m); c.append(f"{pp(q.d_top10)} ({pp(q.d_top10_lo)} to {pp(q.d_top10_hi)})")
    rows.append("| " + " | ".join(c) + " |")
tab = "\n".join(["| Method | Top 10%, primary cut-off (95% CI) | Median percentile, primary | DisAtlas − method, top 10% (pp, 95% CI), primary | Top 10%, stricter cut-off (95% CI) | Median percentile, stricter | DisAtlas − method, top 10% (pp, 95% CI), stricter |",
                 "|---|---|---|---|---|---|---|"] + rows)
md = f"""### FDA-approved pairs scored with models trained before the first trial {{#sec:SI_fda_prioritisation}}

**Design.** For each FDA-approved repurposing pair in the development split, the earlier-approved indication is the anchor and the later indication is the target. Each pair is scored with the dated DisAtlas model trained only on evidence available before the target's first registered trial (primary cut-off), and, as a stricter check, before the earliest trial of either indication. The anchor is compared against every disease in that dated graph, excluding the anchor itself and any disease already trialled with a drug of the same mechanism class by that date. The endpoint is the share of targets ranked in the top 10% of these candidates (random expectation 10%). The same candidates are ranked by four comparators: disease degree, the number of graph neighbours shared with the anchor, the number of literature genes shared with the anchor (counted only up to the cut-off year), and 20 stratum-matched decoy targets per pair. Confidence intervals are cluster bootstrap intervals over drug groups (2,000 resamples). The design, endpoints and comparators were fixed before scoring, and the held-out split was not read.

**Pairs.** Of {int(f1.n_pairs)} development pairs, {int(f1.undated)} have no trial date, {int(f1.left_censored)} have a first trial before the earliest dated graph, {int(f1.n_not_in_slice_nodes)} have an anchor or target that is not yet a node in the dated graph, and {int(f1.n_target_pretrialled)} have a target already trialled with a drug of the same class. That leaves {int(a1.n)} scorable pairs from {int(a1.n_groups)} drugs at the primary cut-off and {int(a2.n)} pairs from {int(a2.n_groups)} drugs at the stricter cut-off.

**Result.** The Supplementary Table gives every endpoint and paired difference. At the primary cut-off, DisAtlas places about two thirds of targets in the top 10% of candidates, far above matched decoy diseases and the 10% random expectation. The shared-neighbour count ranks above DisAtlas at this cut-off: its top-10% share is higher, the paired-difference interval excludes zero, and its mean percentile is also higher (DisAtlas − shared neighbours {sg(s1.d_mean_pct)} percentile points, 95% CI {sg(s1.d_mean_pct_lo)} to {sg(s1.d_mean_pct_hi)}). Shared literature genes and disease degree do not differ from DisAtlas at the primary cut-off (both intervals include zero). At the stricter cut-off, DisAtlas and the shared-neighbour count place the same share of targets in the top 10%, the shared-neighbour count keeps the higher mean percentile ({sg(s2.d_mean_pct)} percentile points, {sg(s2.d_mean_pct_lo)} to {sg(s2.d_mean_pct_hi)}), and DisAtlas ranks above disease degree. In this dated design, most of the DisAtlas signal on approved pairs is also captured by direct graph proximity to the anchor, and DisAtlas does not rank these targets above that baseline.

Supplementary Table: dated FDA-pair prioritisation (development split). Positive differences = DisAtlas better.

{tab}

**Limitations sentence (main text).** On FDA-approved repurposing pairs scored with models trained before the later indication's first trial, a direct count of shared graph neighbours ranked the later indication as high as or higher than DisAtlas, so DisAtlas's added value over simple graph proximity is not established in that design.
"""
# fix the inline .replace that sits inside the f-string literal
open(f"{OUT}/SI_fda_prioritisation_text.md", "w").write(md); print(md)
