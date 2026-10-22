"""Main text, Methods "Comparator systems" -- Table tbl:protocols.

Descriptive cells are fixed protocol text; every count in the table (runs per LLM system,
Phylo seeds, unscored / unmapped pairs) is computed from the stored per-run comparator
scores on the 4,364 IxIDN-Matched pairs (tables/inputs/comparator_run_scores.tsv; LLM
scores are the saved outputs of the runs used in the paper, no model is called here).
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import INPUTS, write_and_check  # noqa: E402

runs = pd.read_csv(INPUTS / "comparator_run_scores.tsv", sep="\t")
assert len(runs) == 4364
systems = {}
for col in runs.columns:
    if "|" in col:
        sysname, run = col.split("|", 1)
        systems.setdefault(sysname, []).append(run)

n_runs = {k: len(v) for k, v in systems.items()}
llm = {"Claude-Sci Opus 5.5": "Claude Sci (Opus 5.5, {} runs)",
       "Claude-Chat Opus 5.5": "Claude Chat (Opus 5.5, {} runs)",
       "Gemini Pro 3.1": "Gemini Pro 3.1 ({} runs)"}
llm_examples = ", ".join(t.format(n_runs[k]) for k, t in llm.items())

phylo_cols = [f"Phylo|{r}" for r in systems["Phylo"]]
seeds = sorted(int(re.search(r"seed(\d+)", r).group(1)) for r in systems["Phylo"])
assert seeds == list(range(seeds[0], seeds[-1] + 1))
phylo_unscored = int(runs[phylo_cols].isna().all(axis=1).sum())
assert phylo_unscored == int(runs[phylo_cols].isna().any(axis=1).sum())
txgnn_unmapped = int(runs["TxGNN|DistMult"].isna().sum())

rows = [
    ["1", "DisAtlas-UPNA-DDAL", "DisAtlas (unseen diseases)",
     "Trained on MechDist with every benchmark disease withheld from embedding and decoder "
     "training; evaluated on the blinded IxIDN-Matched pairs"],
    ["2", "LLM (parametric only)", llm_examples,
     "Structured prompting; DDA score from parametric knowledge only; no tool calls"],
    ["3", "Graph-embedding baseline (not retrained on IxIDN-Matched)",
     f"Phylo (Phylo-DeepWalk, mean of seeds {seeds[0]}-{seeds[-1]})",
     f"Mean score over {['zero','one','two','three','four','five','six'][len(seeds)]} seeds; "
     f"{phylo_unscored} unscored pairs ranked last"],
    ["4", "Knowledge-graph foundation model (published checkpoint)", "TxGNN [@TxGNN]",
     "Published checkpoint, DistMult decoder; not retrained on MechDist or on IxIDN-Matched; "
     f"{txgnn_unmapped} unmapped pairs ranked last"],
]
df = pd.DataFrame(rows, columns=["Protocol", "System type", "Examples", "Evaluation approach"])

EXPECTED = [
    ["1", "DisAtlas-UPNA-DDAL", "DisAtlas (unseen diseases)",
     "Trained on MechDist with every benchmark disease withheld from embedding and decoder "
     "training; evaluated on the blinded IxIDN-Matched pairs"],
    ["2", "LLM (parametric only)",
     "Claude Sci (Opus 5.5, 5 runs), Claude Chat (Opus 5.5, 5 runs), Gemini Pro 3.1 (3 runs)",
     "Structured prompting; DDA score from parametric knowledge only; no tool calls"],
    ["3", "Graph-embedding baseline (not retrained on IxIDN-Matched)",
     "Phylo (Phylo-DeepWalk, mean of seeds 42-46)",
     "Mean score over five seeds; 291 unscored pairs ranked last"],
    ["4", "Knowledge-graph foundation model (published checkpoint)", "TxGNN [@TxGNN]",
     "Published checkpoint, DistMult decoder; not retrained on MechDist or on IxIDN-Matched; "
     "213 unmapped pairs ranked last"],
]

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_protocols", EXPECTED,
                    "Summary of cross-evidence generalization protocols applied across system categories.")
