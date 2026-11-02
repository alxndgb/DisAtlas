"""IxIDN-Matched from the stored QC verdicts and the frozen matching table (Methods, IxIDN-Matched benchmark).

Positives are the QC candidates judged STRONG or PLAUSIBLE; each cell of the seven matching keys keeps
min(n_pos, n_neg) pairs of each class, sampled with seed 0.

    python -m disatlas.benchmarks.ixidn_matched --out DIR [--verdicts PATH]

--verdicts defaults to the stored verdicts data/raw/benchmarks/ixidn_matched/llm_curation_verdicts.tsv (WARM);
pass the output of disatlas.benchmarks.matched_qc (COLD QC) to build from a rerun. --check asserts that the
build equals the shipped ixidn_matched_pairs.tsv and data/derived/fig3_ixidn/qc_verdict_summary.tsv (stored
verdicts) or reports the agreement with them (rerun verdicts).
"""
import argparse
import os

import numpy as np
import pandas as pd

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
QC = os.path.join(ROOT, "data", "raw", "benchmarks", "ixidn_matched_qc")
BM = os.path.join(ROOT, "data", "raw", "benchmarks", "ixidn_matched")
SUMMARY = os.path.join(ROOT, "data", "derived", "fig3_ixidn", "qc_verdict_summary.tsv")
KEEP = ("STRONG", "PLAUSIBLE")
REASON = {"STRONG": "shared-drug evidence judged strong", "PLAUSIBLE": "shared-drug evidence judged plausible",
          "WEAK": "shared-drug evidence judged weak", "SUSPECT": "trial link judged a likely artefact",
          "NOT_IN_GRAPH": "no shared trial drug in the graph (not sent to the language model)"}


def cell_match(cell, y, rng):
    k = np.zeros(len(y), bool)
    for cv in np.unique(cell):
        ip, inn = np.where((cell == cv) & (y == 1))[0], np.where((cell == cv) & (y == 0))[0]
        m = min(len(ip), len(inn))
        if m:
            k[rng.choice(ip, m, replace=False)] = True
            k[rng.choice(inn, m, replace=False)] = True
    return k


def summary(V):
    rows = [dict(verdict=v, decision="kept" if v in KEEP else "removed", reason=REASON[v],
                 pairs=int((V.verdict == v).sum())) for v in REASON]
    rows += [dict(verdict="all", decision=d, reason="", pairs=int(V.verdict.isin(KEEP).eq(d == "kept").sum()))
             for d in ("kept", "removed")]
    r = V[V.rejudged.astype(str).str.lower().eq("true")]
    if len(r):
        agree = int((r.verdict.isin(KEEP) == r.first_pass_verdict.isin(KEEP)).sum())
        rows.append(dict(verdict="repeat check", decision="same keep/drop", reason=f"{agree} of {len(r)}",
                         pairs=len(r), share=round(agree / len(r), 3)))
    return pd.DataFrame(rows)


def build(verdicts, qc_dir=QC):
    V = pd.read_csv(verdicts, sep="\t", keep_default_na=False)
    U = pd.read_parquet(os.path.join(qc_dir, "matching_pairs.parquet"))
    keep = set(map(frozenset, zip(V.disease_a[V.verdict.isin(KEEP)], V.disease_b[V.verdict.isin(KEEP)])))
    pos = np.array([frozenset(k) in keep for k in zip(U.disease_a, U.disease_b)])
    sub = U[(U.y.values == 0) | pos].reset_index(drop=True)
    M = sub[cell_match(sub.cell.values, sub.y.values, np.random.default_rng(0))]
    if int(M.y.sum()) < len(keep):  # cannot happen with the stored verdicts (test_ixidn_matched_qc.py)
        print(f"{len(keep) - int(M.y.sum())} kept positives found no matched negative")
    sp = pd.read_csv(os.path.join(ROOT, "data", "derived", "fig3_ixidn", "disease_split.tsv"), sep="\t")
    sa, sb = M.disease_a.map(dict(zip(sp.mondo_id, sp.split))), M.disease_b.map(dict(zip(sp.mondo_id, sp.split)))
    M = M.assign(y=M.y.astype(int), split=sa.where(sa == sb, "cross"))  # test = both diseases in the test half
    L = pd.read_csv(os.path.join(ROOT, "data", "raw", "benchmarks", "disease_labels.tsv"), sep="\t")
    rare = dict(zip(L.disease_id, L.is_rare.astype(bool)))
    n = M.disease_a.map(rare).astype(int) + M.disease_b.map(rare).astype(int)
    M["pair_class"] = n.map({2: "rare-rare", 1: "rare-common", 0: "common-common"})
    return M[["pair_id", "disease_a", "disease_b", "y", "split", "pair_class"]].reset_index(drop=True), summary(V)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", required=True)
    ap.add_argument("--verdicts", default=os.path.join(BM, "llm_curation_verdicts.tsv"))
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    M, S = build(a.verdicts)
    M.to_csv(os.path.join(a.out, "ixidn_matched_pairs.tsv"), sep="\t", index=False)
    S.to_csv(os.path.join(a.out, "qc_verdict_summary.tsv"), sep="\t", index=False)
    print(f"{len(M)} pairs ({int(M.y.sum())} positive), {int((M.split == 'test').sum())} test")
    if a.check:
        check(M, S, stored=os.path.samefile(a.verdicts, os.path.join(BM, "llm_curation_verdicts.tsv")))


def check(M, S, stored):
    ref = pd.read_csv(os.path.join(BM, "ixidn_matched_pairs.tsv"), sep="\t")
    if stored:
        pd.testing.assert_frame_equal(M.sort_values("pair_id").reset_index(drop=True),
                                      ref.sort_values("pair_id").reset_index(drop=True), check_dtype=False)
        pd.testing.assert_frame_equal(S.fillna("").astype(str), pd.read_csv(SUMMARY, sep="\t", dtype=str,
                                                                            keep_default_na=False))
        print("IxIDN-Matched and the QC summary equal the shipped files")
        return
    a, b = set(M.pair_id[M.y == 1].dropna()), set(ref.pair_id[ref.y == 1])
    print(f"rerun QC: {len(a & b)} of {len(b)} stored positives kept, {int(M.y.sum()) - len(a & b)} new positives")


if __name__ == "__main__":
    main()
