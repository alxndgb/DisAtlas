"""Appendix figure #fig:SI_ixidn_matched_funnel (construction of the IxIDN-Matched benchmark and its shortcut audit).

(a) Disease pairs left after each construction step, for IxIDN- and IxIDN+ (block width proportional to sqrt(pairs)).
(b) Shortest path between the two diseases of a pair in the clinical-trial co-investigation graph (IxIDN+ edges),
    for all 950,131 disease pairs and for the IxIDN-Matched pairs, log scale.

Every count is recomputed from shipped files:
    data/00_raw/frozen/pinned_train_strip_pairs.parquet       IxIDN+ edges (1 hop) and IxIDN- pairs >= 3 hops apart
    data/00_raw/frozen/ixidn_eval_pairs.parquet               IxIDN- pairs that pass the co-investigation veto
    data/00_raw/frozen/ixidn_eval_pair_universe_frozen.tsv.gz evaluation pair universe (985 diseases)
    data/00_raw/frozen/ixidn_matched_universe_features.parquet               universe pairs kept for matching
    data/00_raw/frozen/ixidn_matched_positive_candidates.parquet             positives with a negative in their 7-key cell
    data/00_raw/llm_annotations/ixidn_judge_verdicts.parquet  stored LLM mechanism-QC verdicts
    data/08S_eval/ixidn4364_scores.tsv                        the 4,364 IxIDN-Matched pairs
Writes source_data/ixidn_matched_funnel_counts.tsv and source_data/ixidn_matched_hop_distribution.tsv.

Note on identifiers: the trial-graph pair files write Orphanet diseases as ORPHA:, the universe file as ORPHANET:.
The benchmark was built with these strings as they are, so universe negatives involving an Orphanet disease were
not found in the veto-passing set and were not carried into the matching step (155,375 instead of 158,813). The
script keeps the strings unchanged so that it reproduces the benchmark that was evaluated, and prints both counts.

    python figures/SI/si_ixidn_matched_funnel.py
"""
import os
import sys
from collections import Counter, deque

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Patch, Polygon
from matplotlib.ticker import FixedLocator

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import REPO_ROOT, SOURCE_DATA, apply_compact_style, compact_panel_letter as panel_letter, save

RAW = os.path.join(REPO_ROOT, "data/00_raw")
rd = lambda p: pd.read_parquet(os.path.join(RAW, p))
P = rd("frozen/pinned_train_strip_pairs.parquet")
EV = rd("frozen/ixidn_eval_pairs.parquet")
U = pd.read_csv(os.path.join(RAW, "frozen/ixidn_eval_pair_universe_frozen.tsv.gz"), sep="\t")
F = rd("frozen/ixidn_matched_universe_features.parquet")
Cd = rd("frozen/ixidn_matched_positive_candidates.parquet")
V = rd("llm_annotations/ixidn_judge_verdicts.parquet")
M = pd.read_csv(os.path.join(REPO_ROOT, "data/08S_eval/ixidn4364_scores.tsv"), sep="\t")

key = lambda a, b: (a, b) if a < b else (b, a)
keys = lambda d: set(map(key, d.a, d.b))
pos, negp, negv = P[P["class"] == "pos"], P[P["class"] == "neg"], EV[EV["class"] == "neg"]
Kv, Kpos = keys(negv), keys(pos)
assert Kpos == keys(EV[EV["class"] == "pos"]) and Kv <= keys(negp)

# ---- trial co-investigation graph: BFS hop distance between all disease pairs
nodes = sorted(set(P.a) | set(P.b)); N = len(nodes); ix = {n: i for i, n in enumerate(nodes)}
adj = [[] for _ in nodes]
for a, b in zip(pos.a, pos.b):
    adj[ix[a]].append(ix[b]); adj[ix[b]].append(ix[a])
D = np.full((N, N), -1, dtype=np.int16)
for s in range(N):
    D[s, s] = 0; q = deque([s])
    while q:
        u = q.popleft()
        for w in adj[u]:
            if D[s, w] < 0:
                D[s, w] = D[s, u] + 1; q.append(w)
hop = lambda d: D[[ix[x] for x in d.a], [ix[x] for x in d.b]]

# ---- universe, matching candidates, QC, matched pairs
Uneg = U[U.label == 0]
in_veto = np.array([key(a, b) in Kv for a, b in zip(Uneg.a, Uneg.b)])
Fneg, Fpos = F[F.label == 0], F[F.label == 1]
assert keys(Fneg) == keys(Uneg[in_veto]) and len(Fpos) == int((U.label == 1).sum())
orph = lambda s: s.str.replace("^ORPHA:", "ORPHANET:", regex=True)
Kv_same_prefix = set(map(key, orph(negv.a), orph(negv.b)))
n_same_prefix = sum(key(a, b) in Kv_same_prefix for a, b in zip(Uneg.a, Uneg.b))
qc = V.verdict.value_counts().to_dict()
keep = V[V.verdict.isin(["STRONG", "PLAUSIBLE"])]
assert len(V) == len(Cd) and keys(V.rename(columns={"disease_a": "a", "disease_b": "b"})) == keys(Cd)
Mpos, Mneg = M[M.y == 1], M[M.y == 0]
assert len(Mpos) == len(Mneg) == len(keep) and keys(Mneg) <= Kv
n_dis = len(set(M.a) | set(M.b))

counts = pd.DataFrame([
    ("all", "Disease pairs among trial-graph diseases", N * (N - 1) // 2),
    ("IxIDN+", "Co-investigated (shared trial drug, 1 hop)", len(pos)),
    ("IxIDN-", ">=3 hops apart", len(negp)),
    ("IxIDN-", "Co-investigation veto", len(negv)),
    ("IxIDN+", "Both diseases in the 985-disease universe", len(Fpos)),
    ("IxIDN-", "Both diseases in the 985-disease universe", len(Fneg)),
    ("IxIDN+", "Has a negative in its 7-key shortcut cell", len(Cd)),
    ("IxIDN+", "LLM mechanism QC: STRONG or PLAUSIBLE", len(keep)),
    ("IxIDN+", "IxIDN-Matched", len(Mpos)),
    ("IxIDN-", "IxIDN-Matched (1:1 exact match to a QC-passed positive)", len(Mneg)),
], columns=["class", "step", "pairs"])
counts.to_csv(os.path.join(SOURCE_DATA, "ixidn_matched_funnel_counts.tsv"), sep="\t", index=False)


def hist(h):
    c = Counter(np.where(h < 0, 99, h).tolist()); return c


H = {"all_pairs": hist(D[np.triu_indices(N, 1)]), "ge3_hops": hist(hop(negp)), "veto_passing": hist(hop(negv)),
     "universe_pool": hist(hop(Fneg)), "matched_neg": hist(hop(Mneg)), "matched_pos": hist(hop(Mpos))}
hops_k = sorted(set().union(*H.values()))
hopdf = pd.DataFrame([{"hops": "no path" if k == 99 else k, **{s: H[s].get(k, 0) for s in H}} for k in hops_k])
hopdf.to_csv(os.path.join(SOURCE_DATA, "ixidn_matched_hop_distribution.tsv"), sep="\t", index=False)

# ---- appendix / caption values
c = dict(zip(zip(counts["class"], counts.step), counts.pairs))
assert (N, N * (N - 1) // 2) == (1379, 950131)
assert [len(pos), len(negp), len(negv), len(Fpos), len(Fneg), len(Cd), len(keep), len(Mpos), len(Mneg), n_dis] == \
    [43112, 413856, 335557, 27415, 155375, 4302, 2182, 2182, 2182, 938]
assert qc == {"WEAK": 1354, "STRONG": 1285, "PLAUSIBLE": 897, "SUSPECT": 709, "NOT_IN_GRAPH": 57}
assert {k: v for k, v in H["matched_neg"].items()} == {3: 2100, 4: 80, 5: 2} and dict(H["matched_pos"]) == {1: 2182}
print(f"OK: funnel and hop counts reproduced; universe negatives in the veto set with ORPHA:/ORPHANET: unified: "
      f"{n_same_prefix:,} (benchmark used {len(Fneg):,})")

# ---------------------------------------------------------------- figure
apply_compact_style()
NEG, POS, NEGL, POSL = "#C0392B", "#2E6DB4", "#EBB3AC", "#A9C4E6"
negS = [(c[("all", "Disease pairs among trial-graph diseases")], "All disease pairs\n(1,379 trial-graph diseases)"),
        (len(negp), "≥3 hops apart"), (len(negv), "Co-investigation veto"),
        (len(Fneg), "Both diseases in the\n985-disease universe"), (len(Mneg), "1:1 exact match to a\nQC-passed positive")]
posS = [(len(pos), "Co-investigated\n(shared trial drug, 1 hop)"), (len(Fpos), "Both diseases in the\n985-disease universe"),
        (len(Cd), "Has a negative in its\n7-key shortcut cell"), (len(keep), "LLM mechanism QC:\nSTRONG or PLAUSIBLE"),
        (len(Mpos), "1:1 exact match")]
fig, (a, b) = plt.subplots(2, 1, figsize=(7.2, 7.6), gridspec_kw=dict(height_ratios=[1.25, 1], hspace=0.32))
W = lambda n: 2.3 * np.sqrt(n / 950131) + 0.09


def funnel(st, cx, col, light, side, lx):
    for i, (n, lab) in enumerate(st):
        top = W(n); bot = W(st[i + 1][0]) if i + 1 < len(st) else W(n)
        yt = -i; yb = yt - 0.88; ym = (yt + yb) / 2
        a.add_patch(Polygon([(cx - top / 2, yt), (cx + top / 2, yt), (cx + bot / 2, yb), (cx - bot / 2, yb)],
                            closed=True, fc=light if i < len(st) - 1 else col, ec="white", lw=1))
        if top > 0.7:
            a.text(cx, ym, f"{n:,}", ha="center", va="center", fontsize=8, fontweight="bold")
        else:
            a.text(cx - side * (max(top, bot) / 2 + 0.08), ym, f"{n:,}", ha="right" if side > 0 else "left",
                   va="center", fontsize=8, fontweight="bold", color=col)
        a.text(lx, ym, lab, ha="right" if side > 0 else "left", va="center", fontsize=7, color="#333")


# panel b
hops = ["no path", "1", "2", "3", "4", "5", "6"]
order = [99, 1, 2, 3, 4, 5, 6]
allp = [H["all_pairs"].get(k, 0) for k in order]
matched = [H["matched_pos"].get(k, 0) + H["matched_neg"].get(k, 0) for k in order]
cl = ["#C4C2D1", POSL, "#9A9A9A", NEGL, NEGL, NEGL, NEGL]; cd = [None, POS, None, NEG, NEG, NEG, NEG]
x = np.arange(7)
b.bar(x, allp, 0.8, color=cl)
b.bar(x, [m if m else np.nan for m in matched], 0.45, color=[k or "none" for k in cd])
for i, v in enumerate(allp): b.text(i, v * 1.3, f"{v:,}", ha="center", fontsize=6)
for i, m in enumerate(matched):
    if m > 50: b.text(i, m / 1.6, f"{m:,}", ha="center", va="top", fontsize=6, color="white", fontweight="bold")
    elif m: b.text(i + 0.3, m * 1.2, f"{m:,}", ha="left", va="bottom", fontsize=6, color=NEG, fontweight="bold")
b.axvline(2.5, ls="--", color="#555", lw=.8); b.text(2.57, 2.2e7, "IxIDN⁻ sampled from ≥3 hops", fontsize=6, va="top", color=NEG)
b.set_yscale("log"); b.set_ylim(1, 4e7); b.set_xticks(x); b.set_xticklabels(hops)
b.yaxis.set_major_locator(FixedLocator([1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7]))
b.set_yticklabels(["1", "10", "100", "1k", "10k", "100k", "1M", "10M"])
b.set_xlabel("Shortest path between the two diseases in the trial co-investigation graph (hops)")
b.set_ylabel("Disease pairs (log scale)")
b.set_title("Hop distance: matched positives are 1 hop, matched negatives ≥3 hops", loc="left")
b.legend(handles=[Patch(fc="#9A9A9A", label="All disease pairs (light shades)"), Patch(fc=NEG, label="IxIDN-Matched IxIDN⁻"),
                  Patch(fc=POS, label="IxIDN-Matched IxIDN⁺")], frameon=False, fontsize=6, loc="upper right")
panel_letter(b, "b")

# panel a
funnel(negS, -1.55, NEG, NEGL, 1, -2.85); funnel(posS, 1.3, POS, POSL, -1, 2.75)
a.text(-1.55, 0.3, "IxIDN⁻ (negatives)", ha="center", fontsize=8, color=NEG, fontweight="bold")
a.text(1.3, 0.3, "IxIDN⁺ (positives)", ha="center", fontsize=8, color=POS, fontweight="bold")
a.add_patch(FancyBboxPatch((-3.65, -5.95), 6.8, 0.62, boxstyle="round,pad=0.03", fc="#F3F3F3", ec="#888", lw=.8))
a.text(-0.25, -5.64, f"IxIDN-Matched: {len(M):,} pairs ({len(Mpos):,} IxIDN⁺ + {len(Mneg):,} IxIDN⁻), {n_dis} diseases\n"
       "out-of-distribution transfer test", ha="center", va="center", fontsize=7)
for cx, tx in ((-1.55, -1.1), (1.3, 0.7)):
    a.annotate("", xy=(tx, -5.3), xytext=(cx, -4.9), arrowprops=dict(arrowstyle="-|>", color="#888", lw=.8))
a.set_xlim(-5.0, 5.4); a.set_ylim(-6.05, 0.55); a.axis("off")
a.set_title("How IxIDN-Matched is built: pairs left after each step (width ∝ √pairs)", loc="left")
panel_letter(a, "a")
fig.subplots_adjust(hspace=0.06)
pa, pb = a.get_position(), b.get_position(); a.set_position([pa.x0, pb.y1 + 0.06, pa.width, pa.height])
print(save(fig, "ixidn_matched_funnel", bbox_inches="tight"))
