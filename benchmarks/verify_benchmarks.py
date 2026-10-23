"""Rebuild the comparator benchmark numbers from stored outputs and assert they equal the manuscript.

    python benchmarks/verify_benchmarks.py

Default stored-output path only: no network access, no model calls, no key.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
for p in (HERE, os.path.join(HERE, "llm"), os.path.join(HERE, "txgnn")):
    sys.path.insert(0, p)
import score_llm  # noqa: E402
import score_txgnn  # noqa: E402

# Manuscript values at reported precision.
AUROC = {  # Results (cross-system comparison) and Supplementary "LLM subgroups"
    ("all pairs", "DisAtlas (unseen in learning)"): 0.742, ("all pairs", "DisAtlas (full)"): 0.765,
    ("all pairs", "Claude Chat (Opus 5.5)"): 0.692, ("all pairs", "Claude Sci (Opus 5.5)"): 0.685,
    ("all pairs", "Gemini Pro 3.1"): 0.497,
    ("held-out test split", "DisAtlas (unseen in learning)"): 0.754,
    ("held-out test split", "Claude Chat (Opus 5.5)"): 0.696, ("held-out test split", "Claude Sci (Opus 5.5)"): 0.702,
    ("held-out test split", "Gemini Pro 3.1"): 0.500,
    ("non-oncology", "DisAtlas (unseen in learning)"): 0.745,
    ("non-oncology", "Claude Chat (Opus 5.5)"): 0.671, ("non-oncology", "Claude Sci (Opus 5.5)"): 0.664,
    ("non-oncology", "Gemini Pro 3.1"): 0.428,
    ("oncology", "DisAtlas (unseen in learning)"): 0.718, ("oncology", "Claude Sci (Opus 5.5)"): 0.724,
    ("oncology", "Claude Chat (Opus 5.5)"): 0.701, ("oncology", "Gemini Pro 3.1"): 0.559,
}
N = {"all pairs": 4364, "oncology": 1740, "non-oncology": 2624}
AGREE = {  # Figure 6c: median (min-max) pairwise Spearman
    "DisAtlas (unseen in learning)": (0.973, 0.971, 0.976), "DisAtlas (full)": (0.921, 0.912, 0.928),
    "Claude Sci (Opus 5.5)": (0.645, 0.363, 0.831), "Claude Chat (Opus 5.5)": (0.533, 0.152, 0.807),
    "Gemini Pro 3.1": (0.555, 0.504, 0.873),
}


def main():
    d, A, AGs, MT, J, onc = score_llm.main([])
    A = A.set_index(["subset", "system"])
    for key, v in AUROC.items():
        assert round(A.loc[key, "auroc"], 3) == v, (key, A.loc[key, "auroc"], v)
    for s, n in N.items():
        assert int(A.loc[s].n.iloc[0]) == n, s
    AGs = AGs.set_index("method")
    for m, v in AGREE.items():
        got = tuple(round(AGs.loc[m, c], 3) for c in ("median", "min", "max"))
        assert got == v, (m, got, v)
    # oncology: DisAtlas (unseen) and Claude Sci "statistically indistinguishable"
    assert onc.lo.item() < 0 < onc.hi.item()
    _, _, cov, TM = score_txgnn.main()
    assert int(cov.ranked_last.sum()) == 213
    assert round(TM.set_index("metric").loc["AUROC", "value"], 3) == 0.566
    print("benchmarks: all manuscript values reproduced")


if __name__ == "__main__":
    main()
