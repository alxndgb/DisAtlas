| Protocol | System type | Examples | Evaluation approach |
|---|---|---|---|
| 1 | DisAtlas-UPNA-DDAL | DisAtlas (unseen diseases) | Trained on MechDist with every benchmark disease withheld from embedding and decoder training; evaluated on the blinded IxIDN-Matched pairs |
| 2 | LLM (parametric only) | Claude Sci (Opus 5.5, 5 runs), Claude Chat (Opus 5.5, 5 runs), Gemini Pro 3.1 (3 runs) | Structured prompting; DDA score from parametric knowledge only; no tool calls |
| 3 | Graph-embedding baseline (not retrained on IxIDN-Matched) | Phylo (Phylo-DeepWalk, mean of seeds 42-46) | Mean score over five seeds; 291 unscored pairs ranked last |
| 4 | Knowledge-graph foundation model (published checkpoint) | TxGNN [@TxGNN] | Published checkpoint, DistMult decoder; not retrained on MechDist or on IxIDN-Matched; 213 unmapped pairs ranked last |

Summary of cross-evidence generalization protocols applied across system categories.
