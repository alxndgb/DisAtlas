# DisAtlas release report (decoder-free, cosine-only)

Binding text: the final story markdown (sha256 prefix `db3b7eb113f4`, commit c0b04efc, with the author's Table S13 edit; it differs from the earlier `52e96bb6ce12` only in the wording of lines 683, 692 and 2003 and the four S13 values on line 1737): 7 main figures, 26 Supplementary figures, Table 1 and
18 Supplementary tables. Every DisAtlas pair number is the mean over seeds 4707-4711 of the cosine of the two 32-d embeddings;
the clustering distance is 1 - cosine. The Hadamard-MLP decoder is not shipped.

## 1. Rebuild classes of derived inputs (`release/rebuild_classes.tsv`)

23 items: 15 class **a** (rebuilt from released inputs, equal at stored precision) and 8 class **c** (rebuilt but not
bit-identical, or not rebuilt; shipped values kept). No class b remains. Each row names the commit and test
(18 `tests/test_rebuild_*.py`). Class c items:

| item | reason |
|---|---|
| Figure 7 (not rebuilt) | per-event comparison-set sizes and per-candidate lists are not released |
| Figure S16 c | exact except three within-cluster pair-share estimates (argsort tie order) |
| Figure S18 | float32 vs float64 near-tie at MONDO:0013751; equal at the 3 d.p. printed |
| Table S8 | DisAtlas and node-disjoint rows: released 5-seed cosines vs internal per-pair scores (AUROC +0.002 to +0.013); all other rows exact |
| Table S10 (two rows) | the authors' disease-gene table includes 610 pairs that are not redistributed (604 supported by OMIM alone, 6 with no public source) |
| Figure S3 / S12 text, HPO similarity | restricted HPO evidence |
| Figure S13 legacy tables | decoder-era, removed (section 3) |

## 2. For the author: printed values that released inputs do not reproduce

1. **Table S13 / SI top-10 precision** (author decision 2026-11-15; `release/text_edits_for_author.tsv`): settled. The
   shipped table is the rebuild from released inputs (gene-anchored models, cosine readout, the run behind every other
   event-time number); the comparator columns are unchanged. The four DisAtlas means were printed from an earlier scoring
   run; the author has updated the text to the released values, and the parity checks are bound to them (the S13
   `manuscript_mismatch` list is empty). Stated orderings still hold.

   | registry id | location | earlier printed | released (shipped, now printed) | status |
   |---|---|---|---|---|
   | SI05372 | prec_s4 test, DisAtlas mean | 0.57% | 0.53% | applied by author (2026-11-15) |
   | SI05375 | prec_s4 all, DisAtlas mean | 0.37% | 0.35% | applied by author (2026-11-15) |
   | SI05376 | prec_s3 test, DisAtlas mean | 9.91% | 10.21% | applied by author (2026-11-15) |
   | SI05379 | prec_s3 all, DisAtlas mean | 5.50% | 5.79% | applied by author (2026-11-15) |

2. **SI "Cluster themes" subsection (Table S10)**: the two-sided Mann-Whitney p for gene sharing ("7.8-fold against
   17.2-fold; p = 0.032") becomes 0.052 with the released inputs, because the gene/GO pair-overlap columns are class c.
   Reproducing the printed value needs the 610 non-redistributed pairs: removing them alone (and nothing else) moves p from
   0.032 to 0.052 (controlled substitution of each input difference; data/README.md, "Table S10 (cluster themes)
   and the disease-gene table").
3. **Headline test AUROC**: the shipped `headline_test_auroc.tsv` held a stale 0.859 [0.803-0.905]. The SI numbers script now
   reads `fig4_values.json`, which holds the printed cosine values, and the stale file is removed. No asserted parity key
   changed.
4. **Sibio parity rows still to decide**: `release/sibio_review_for_lead.tsv` lists 68 non-bound rows (64 unchecked kept,
   1 moved to unchecked, 3 pinned unprinted), each with an empty `lead_decision` column (keep / bind / retire).
5. **Licences to confirm** (authors' internal licence table, not in the repository; `release/LICENSE_PROPOSAL.md`): set
   by the author on 2026-11-15 and still marked "author to confirm" (not verified on the sources): HPO = HPO licence
   (attribution and version; ontology content not altered, propagated edges are derived data); Zenodo `event_time/` =
   CC BY-SA 3.0 (ChEMBL share-alike); training graph, BIKG HABCD subset and restricted-source salted tokens = CC BY 4.0
   for the derived salted edges only, with data-owner clearance pending. Still open: Mondo and GO release versions, the
   remaining public rows, stored language-model outputs. The code licence is MIT (`LICENSE`, author decision 2026-11-15);
   copyright holder: The DisAtlas authors.

## 3. Removals (`release/MINIMISATION.tsv`, files in the trash, nothing deleted)

The removed files are kept in `/scratch/kwfh644/DisAtlas_archive_1115/DisAtlas_v32_trash/step4_prune_20261115/`.

Files removed in this round. A file-open audit of the full pipeline and test suite, plus a grep of code, notebooks, tests
and registries, found no reader for any of them:

- `data/raw/benchmarks/ordo/ixidn_matched_ontology_distance.tsv` (66,398 B) and `ordo_distance_summary.tsv` (167 B). The
  ontology baseline reads `ixidn4364_ontology_distance.tsv`; the staging block that wrote them is removed.
- `data/derived/si_rare_common/headline_test_auroc.tsv` (80 B, stale; see 2.3).
- `data/derived/uncertainty/ci_pair_vs_cluster.tsv` (6,062 B) and `headline_core_values.tsv` (4,085 B). These are
  decoder-era tables; notebook 03 now reads `headline_cosine_ci.tsv`.

The six earlier row and column reductions are listed in the same file.

## 4. Zenodo

Not yet deposited: no Zenodo record exists. The manifest and the folder names are final once the author approves them.

`release/ZENODO_MANIFEST.tsv` covers 38 files, 949,726,744 B, giving size, sha256, the checksum lists that cover each file, and
the licence. It is generated by `scripts/maintainer/zenodo_manifest.py` (`--check` returns 0). `data/ZENODO_SHA256SUMS.txt`
now covers all 38 files: 24 entries were added (graph, ontology_term_term_distance [named specific_channels until 2026-11-15], two model_outputs docs). The README check command
was fixed.

## 4a. Trial identifiers and dates (author decision 2026-11-15)

Trial identifiers (NCT ids) and trial dates in `tables/SI`, `figures/SI` and the derived tables they read are released as
shipped. The earlier 'no trial identifiers, no dates' rule is withdrawn. No test or check enforced it; it was stated in the
local agent instructions only. No data file was changed. Notes: `data/README.md` (Licence and blinding) and the
authors' internal licence table (ClinicalTrials.gov row). The code licence is MIT (author decision 2026-11-15; `LICENSE`).

## 4b. Figure 1 and open staging items

Figure 1 is a generated item: `figures/main/fig1_knowledge_transfer.py` draws it from the five tracked files in
`data/derived/fig1_open_targets/` (Open Targets 26.09) with the `figure` subcommand of `src/disatlas/fig1_open_targets.py`
(coverage_table: reproduced, identical bytes; `data/DATA_INDEX.tsv` and `data/README.md` name this script as the reader).
The `tables` subcommand builds those five files from the Open Targets release; its optional first step,
`scripts/fetch_open_targets.py` (public Open Targets GraphQL API and bulk parquet only), is staged from the private repository
(read in full: no private paths, credentials, internal hosts or internal data reads; its data root comes from `--data-root`,
`$OT_DATA_ROOT` or `./open_targets_data`, and the private `fetch_config.json` is not shipped). One edit: `requests` is imported
inside the client, so import and `--help` need only requirements.txt (`fetch` needs `pip install requests`). It is not run by
`run_pipeline.sh`; `tests/test_fetch_open_targets.py` checks import and `--help` with network access blocked.

Copyright holder in `LICENSE`: The DisAtlas authors (author decision; no longer open).

## 5. Gates (fresh clone, fresh venv, SLURM rd-bds)

Scripts are in `/scratch/kwfh644/DisAtlas_r6_work/gates_1115/`. Each run builds a fresh clone, a fresh Python 3.11 venv
from `requirements.txt` (34 packages) and a copy of `data/zenodo`.

- **WARM**, on dfc5bd3 (the last code or data commit; the report commit that follows adds documentation only):
  `run_pipeline.sh` rc 0, 0 changed tracked files, `pytest tests` gives 3,243 passed, 220 skipped (slow tests, run with
  `DISATLAS_SLOW=1`), 0 failed, in 791 s. The Figure S14 test, which had failed once in an earlier run, passed.
- **COLD**, on 7aeb590 (adds `requirements-retrain.txt`, torch 2.9.1+cpu): the graph matches its recorded signature, all 5
  seeds were retrained, `test_cold_exact.py` and `test_cold_tolerance.py` passed, the re-scored test AUROC is 0.8578
  (released model 0.8605; PARITY rule in README), and there are 0 changed tracked files. The only failure was
  `test_data_index` (stale index, the slip in section 8, fixed in dfc5bd3). Commits between 7aeb590 and dfc5bd3 change
  only documentation and the data index. Total 1,358 s.

## 6. Numbering map

`release/NUMBERING_FINAL.tsv` maps every final item to its crossref id and caption excerpt, in the binding text's order.
Relative to the 2026-11-10 audit snapshot (`release/audit_20261110`), old Figures S16 (negative-set realism) and S28
(compositionality) were removed. Final Figures S1-S15 are old S1-S15, and final S16-S26 are old S17-S27. Table numbers
follow the binding text; use the crossref id rather than old table numbers. Table S18 (`SI_trace_order`) was S22 in the
old numbering. Its companion values table (old S23) is not printed in the final SI and remains a README output only.

## 7. Archive

The work folders `DisAtlas_r4_plan`, `DisAtlas_r5_fresh`, `DisAtlas_r5_work` and `DisAtlas_v32_trash` were moved to
`/scratch/kwfh644/DisAtlas_archive_1115`. File counts were equal before and after each move. Every git bundle (7, all pass
`git bundle verify`) and the three earlier reports were first copied to `/projects/rd/us/AI/kwfh644/DisAtlas_history`
(checksums in its `SHA256SUMS.txt`). An empty `DisAtlas_v32_trash` was recreated for later removals.
`DisAtlas_v32_private` (restricted inputs, blinding salt) was left in place.

## 8. Rule slips

- Once, a one-line Python decompress ran on the login node. Later uploads used `gzip -d`.
- An empty `data/raw/benchmarks/ordo/` folder was removed with `rmdir` after its two files were moved to the trash.
- `data/DATA_INDEX.tsv` was regenerated before `zenodo_manifest.py` was tracked, so the index missed one reader. The
  fresh-clone warm gate caught it, and it was fixed in dfc5bd3.
- A first draft of the sibio review script would have rewritten the registry to strip newlines. It was made read-only
  before it ran.
- All commits are local on `main`. Nothing was pushed.

## 9. Git history (author decision 2026-11-15)

- By author decision the repository history is one root commit with no parent. Its tree is the 11-commit squashed `main`
  (tip `74343fe`) plus this section, without the internal licence table, and with the edits that remove the references to it.
- The authors' internal licence table is not in the repository: it is untracked and listed in `.gitignore`. The data licences are stated in `data/README.md` (Licence and blinding), the code licence is `LICENSE` (MIT).
  Older commits that contain it exist only in the private bundles below, which are not pushed.
- Rollback points, as verified bundles in `/projects/rd/us/AI/kwfh644/DisAtlas_history/bundles/release_squash_1115/` (checksums in
  `DisAtlas_history/SHA256SUMS.txt`):
  - old `main` = rollback point: 147 commits (145 not on `origin/main`), tip `c3d9f7d`, `DisAtlas_main_c3d9f7d.bundle`;
  - 11-commit squashed `main`, tip `74343fe`, `DisAtlas_main_squashed11_74343fe.bundle`.
- Commit hashes quoted above in this report (for example dfc5bd3, b707cf8) refer to the old history and resolve in the
  first bundle.
- The 11-commit file list `release_squashed_structure.tsv` (kept in the release hand-off folder) is now historical.
- `origin/main` is unchanged: two commits that local `main` never had (6543348, an empty 'Initial commit', and 4dc4e60,
  README only), unrelated to the release history. Author decision 2026-11-15: replace the remote with the single local
  root commit; no merge, no rebase onto them. The old remote tip 4dc4e60 is restorable from
  `DisAtlas_origin_main_4dc4e60.bundle` (same folder, in `SHA256SUMS.txt`) or by fetching that sha.
- Until the force-push, `git status` reports `main` as ahead 1, behind 2 of `origin/main`. This is expected and does not
  mean commits are missing; do not pull or merge, and do not run `git fetch --prune` (it would drop the `origin/main` ref
  that the lease below names).
- Push procedure, pending a separate explicit author go-ahead (not executed): over the repository's deploy key,
  `git push --force-with-lease=main:4dc4e60 origin main`. The lease makes the push fail if the remote has moved from 4dc4e60.
