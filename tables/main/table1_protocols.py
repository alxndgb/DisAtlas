"""Main Table 1: protocols of the cross-system comparison on IxIDN-Matched (Methods, Comparator systems).

The table is descriptive text. Run counts and the number of TxGNN-unmapped pairs are derived from the stored per-run
score file so that the table cannot drift from the data behind Figure 6.
    python tables/main/table1_protocols.py      -> table1_protocols.tsv, table1_protocols.md, table1_protocols_values.tsv
"""
import os, sys
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
from disatlas.parity import write_values

D = os.path.join(ROOT, "data", "derived", "fig6_comparison")
R = pd.read_csv(os.path.join(D, "per_run_scores.tsv"), sep="\t")
ND = pd.read_csv(os.path.join(D, "node_disjoint_per_seed.tsv"), sep="\t")
runs = lambda s: sum(c.split("|")[0] == s for c in R.columns if "|" in c)
n = {"claude_chat_runs": runs("Claude-Chat Opus 5.5"), "gemini_runs": runs("Gemini Pro 3.1"),
     "claude_sci_runs": runs("Claude-Sci Opus 5.5"), "phylo_seeds": runs("Phylo"),
     "node_disjoint_seed_columns": sum("|" in c for c in ND.columns),
     "txgnn_unmapped_pairs": int(R["TxGNN published ckpt|DistMult"].isna().sum()), "pairs": len(R)}

rows = [
 (1, "DisAtlas", "DisAtlas (benchmark diseases held out of training, node-disjoint)",
  "Benchmark diseases withheld at every step of network construction and training, five seeds per fold; "
  "evaluated on the blinded IxIDN-Matched pairs"),
 (2, "Large language model", f"Claude Chat (Opus 5.5, {n['claude_chat_runs']} runs), Gemini Pro 3.1 ({n['gemini_runs']} runs)",
  "Structured prompt with the blinded pair file; mean score over runs"),
 (3, "Agentic AI system", f"Claude Sci (Claude Science, Opus 5.5, {n['claude_sci_runs']} runs); Phylo Biomni",
  "Mean score over runs or seeds"),
 (4, "Knowledge-graph foundation model (published checkpoint)", "TxGNN",
  f"Published checkpoint, DistMult decoder; not retrained on MechDist or on IxIDN-Matched; "
  f"{n['txgnn_unmapped_pairs']} unmapped pairs ranked last"),
]
T = pd.DataFrame(rows, columns=["Protocol", "System type", "Examples", "Evaluation approach"])
T.to_csv(os.path.join(HERE, "table1_protocols.tsv"), sep="\t", index=False)
with open(os.path.join(HERE, "table1_protocols.md"), "w") as fh:
    fh.write("| " + " | ".join(T.columns) + " |\n|" + "---|" * len(T.columns) + "\n")
    for r in T.itertuples(index=False): fh.write("| " + " | ".join(map(str, r)) + " |\n")
    fh.write("\nSummary of cross-evidence generalization protocols applied across system categories.\n")
write_values(os.path.join(HERE, "table1_protocols_values.tsv"), n)
print(n)
