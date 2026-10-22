# Language-model comparators (Protocol 2)

    python benchmarks/llm/score_llm.py            # default: stored outputs only

## Stored outputs used in the paper (`data/`)
All tables are keyed by `blind_id` (the opaque pair identifier the systems saw) and have 4,364 rows with no missing scores.
- `protocol2_upload.tsv.gz` is the file uploaded to every system, containing `blind_id`, the MONDO ids and the disease names. Pairs are shuffled and carry no labels.
- `llm_run_scores.tsv.gz` holds the returned DDA score per repeat run, for Claude Chat (Opus 5.5) and Claude Sci (Opus 5.5) with 5 runs each and Gemini Pro 3.1 with 3 runs. Each run was a fresh session with no tool use.
- `disatlas_unseen_per_seed.tsv.gz` holds the DisAtlas scores with the benchmark diseases unseen in learning, one column per training seed (5 seeds). The full-DisAtlas scores come from `data/08S_eval/`.

## Prompts and models
- `prompts/protocol2_prompt.txt` is the prompt used for every Protocol 2 run.
- `models.json` lists the systems, run counts and the prompt md5.

## Scoring
A system's score is its mean over repeat runs (for DisAtlas, the mean over seeds). Unscored pairs are ranked last.
- Subgroup comparison: disease-clustered bootstrap with B=1,000.

## Optional re-run (`--rerun-llm`)
`python benchmarks/llm/score_llm.py --rerun-llm` (or `DISATLAS_RERUN_LLM=1`) calls `rerun_llm.py`, which re-sends the
shipped prompt and upload to the models with your own key. Set `ANTHROPIC_API_KEY` and/or `GEMINI_API_KEY`, and the
Gemini model ID in `models.json`. A key is needed only for this switch. New outputs go to `rerun/` (gitignored) and
are scored there. Model outputs are not deterministic, and the paper runs used each system's own interface with a
file upload, so re-run numbers will differ from the paper.
