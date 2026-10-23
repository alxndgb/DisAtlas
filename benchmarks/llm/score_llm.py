"""Score the language-model comparator runs on IxIDN-Matched from the stored outputs used in the paper.

Manuscript: Results, cross-system comparison (Protocol 2, Table "protocols"); Figure 6a-d;
Supplementary "LLM subgroups".

    python benchmarks/llm/score_llm.py               # default: stored outputs only, no network, no key
    python benchmarks/llm/score_llm.py --rerun-llm   # opt-in: re-query the models first (see rerun_llm.py)

The default path reads only files under benchmarks/llm/data/ and data/08S_eval/. With --rerun-llm the
prompts are re-sent with your own API key, the new DDA columns are written to benchmarks/llm/rerun/, and
those are scored instead. Re-runs are not deterministic and will not reproduce the paper's numbers exactly.
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
from _common import (REPO_ROOT, load_key, read_tsv, ranked_last, auc, pair_bootstrap_metrics, topk_jaccard,  # noqa: E402
                     run_agreement, disease_cluster_bootstrap, write_tsv)

DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "results")

SYSTEMS = {  # display name -> stored run columns
    "Claude Chat (Opus 5.5)": [f"claude_chat_run{i}" for i in range(1, 6)],
    "Claude Sci (Opus 5.5)": [f"claude_sci_run{i}" for i in range(1, 6)],
    "Gemini Pro 3.1": [f"gemini_pro_3_1_run{i}" for i in range(1, 4)],
}
UNSEEN = "DisAtlas (unseen in learning)"
FULL = "DisAtlas (full)"
SEEDS = ["p_s4707", "p_s11", "p_s202", "p_s3407", "p_s8191"]


def load(run_scores=None):
    k = load_key()
    runs = read_tsv(run_scores or os.path.join(DATA, "llm_run_scores.tsv.gz"))
    uns = read_tsv(os.path.join(DATA, "disatlas_unseen_per_seed.tsv.gz")).rename(columns={s: s + "_u" for s in SEEDS})
    d = k.merge(runs, on="blind_id", how="left").merge(uns, on="blind_id", how="left")
    assert len(d) == 4364 and d[[s + "_u" for s in SEEDS]].notna().all().all()
    return d


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rerun-llm", action="store_true",
                    help="opt-in: re-query the models with your own API key before scoring (non-deterministic)")
    a = ap.parse_args(argv)
    if a.rerun_llm or os.environ.get("DISATLAS_RERUN_LLM") == "1":
        from rerun_llm import rerun_all
        paths = rerun_all()
        d = load(paths["runs"]); out = os.path.join(HERE, "rerun", "results")
    else:
        d = load(); out = OUT
    y = d.y.values.astype(int)

    # ---- system scores: multi-run systems use the mean score across runs; unscored pairs ranked last ----
    raw = {FULL: d.p.values, UNSEEN: d[[s + "_u" for s in SEEDS]].mean(axis=1).values}
    raw.update({nm: d[c].mean(axis=1).values for nm, c in SYSTEMS.items()})
    S = {nm: ranked_last(v) for nm, v in raw.items()}

    # ---- Results / Supplementary subgroups: AUROC by subset ----
    subsets = {"all pairs": np.ones(len(d), bool), "held-out test split": (d.pool == "test").values,
               "oncology": (d.either_is_cancer == 1).values, "non-oncology": (d.either_is_cancer == 0).values}
    rows = [dict(subset=sn, n=int(m.sum()), system=nm, auroc=auc(y[m], S[nm][m]))
            for sn, m in subsets.items() for nm in S]
    A = pd.DataFrame(rows); write_tsv(A, os.path.join(out, "protocol2_auroc_by_subset.tsv"))

    # oncology: DisAtlas (unseen) vs Claude Sci, paired disease-clustered bootstrap
    on = d[subsets["oncology"]].copy(); on["u"] = S[UNSEEN][subsets["oncology"]]
    on["c"] = S["Claude Sci (Opus 5.5)"][subsets["oncology"]]
    pt, bs = disease_cluster_bootstrap(on, ["u", "c"], n_boot=1000)
    dl = bs["u"] - bs["c"]
    onc = pd.DataFrame([dict(subset="oncology", contrast=f"{UNSEEN} minus Claude Sci (Opus 5.5)",
                             delta=pt["u"] - pt["c"], lo=np.quantile(dl, .025), hi=np.quantile(dl, .975))])
    write_tsv(onc, os.path.join(out, "protocol2_oncology_contrast.tsv"))

    # ---- Figure 6c: run-to-run agreement ----
    AG = pd.concat([run_agreement(d, [s + "_u" for s in SEEDS], UNSEEN), run_agreement(d, SEEDS, FULL)]
                   + [run_agreement(d, c, nm) for nm, c in SYSTEMS.items()], ignore_index=True)
    write_tsv(AG, os.path.join(out, "run_agreement.tsv"))
    AGs = AG.groupby("method", sort=False).spearman.agg(median="median", min="min", max="max").reset_index()
    write_tsv(AGs, os.path.join(out, "run_agreement_summary.tsv"))

    # ---- Figure 6b/6d: ranking metrics (pair bootstrap) and top-10% overlap ----
    MT = pair_bootstrap_metrics(y, S); write_tsv(MT, os.path.join(out, "ranking_metrics.tsv"))
    J = topk_jaccard(raw); J.index.name = "method"
    write_tsv(J.reset_index(), os.path.join(out, "top10_jaccard.tsv"))

    print(A.pivot(index="system", columns="subset", values="auroc").round(3).to_string())
    print(AGs.round(3).to_string(index=False))
    print(onc.round(3).to_string(index=False))
    return d, A, AGs, MT, J, onc


if __name__ == "__main__":
    main()
