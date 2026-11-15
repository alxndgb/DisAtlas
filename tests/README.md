# Number-parity tests: source text

The registries `main_text_numbers.tsv` and `si_text_numbers.tsv` bind every number of the final revised manuscript and
SI appendix to a shipped quantity (`main_quantities.py`, `si_quantities.py`, `si_derived_quantities.tsv`).
The binding source is the author's story markdown, which carries the main text and the Supplementary Information and is
always the latest source text; the compiled docx is used only to cross-check figure/table numbering and images.

| file (in the author's manuscript repository) | git commit | sha256 | mtime | size (bytes) |
|---|---|---|---|---|
| `assets/manuscripts/disatlas/disatlas_representation_learning_story.md` (binding source) | `c0b04efc3ed69b09c168cbb036db726e8c3405d4` (2026-11-14 21:59:40 -0500) | `db3b7eb113f41106ad5356318f640124d23fb5b4d10f33b724a99fb84a451749` | 2026-11-14 21:59:29 -0500 | 383,468 |
| `assets/manuscripts/disatlas/output/disatlas_revised_manuscript.docx` (numbering/image cross-check) | build output | `3213851cfbf98f02a39163d9020994f55cb88bbb98f2c4720b186248d87558dd` | 2026-11-14 12:16:42 -0500 | 16,684,730 |
| `revision/response/output/response_to_referees.docx` | build output | `fa65816ca8ac1c1abecd8953a7a79f6c0aec3728abdd82b8937c5486d47fa9cb` | 2026-11-14 11:16:04 -0500 | 1,848,693 |

The registries were bound at commit `2fee522e` (sha256 `52e96bb6ce12`, recorded in the parity specs and review files); commit
`c0b04efc` applies the author's Table S13 edit and changes only lines 683, 692 and 2003 (wording, no numbers) and the four
DisAtlas means on line 1737 (rows SI05372, SI05375, SI05376, SI05379, now `computed`).

Version labels are not stored in either file. The registries were bound to the plain-text extraction of the docx; the
markdown (pandoc markup, citations and cross-references removed) carries the same numbers: every numeric passage of
one is found, digit for digit, in the other. The only one-sided numbers are the YAML metadata of the markdown and the
rendered equation, cross-reference and bibliography numbers of the docx, none of which is a registry value.

Pair numbers are cosine scores: the cosine of the two 32-d embeddings, mean of the five per-seed cosines (seeds
4707-4711); all-pairs percentiles are over the 84,766,710 upper-triangle pairs. Values must match at the printed precision; the only tolerance is the decimal rule of `tolerance.py`
(op `dtol`: |printed - regenerated| <= 0.003 in the printed unit, decimal tokens only, never counts, p-values,
whole percents or odds ratios; each such row gives both values in its note).

Run: `pytest tests/test_main_parity.py tests/test_si_parity.py` (rd-bds; never on a login node).

## Rebuild tests (`test_rebuild_*.py`)

Each `test_rebuild_<item>.py` recomputes a shipped table in `data/derived` from released inputs with
`src/disatlas/rebuild/<item>.py` and compares it with the shipped file (tolerances and any documented exception are stated in
the test docstring). `release/rebuild_classes.tsv` lists the 18 rebuild tests with the item each covers and its rebuild
class (a: equal at stored precision; c: not bit-identical, shipped values kept); the Figure S14 test runs without a flag. Tests that run the full bootstrap or permutation count are skipped unless `DISATLAS_SLOW=1`; the
power-law and K-selection bootstraps also accept the output folder of the SLURM arrays in `scripts/slurm/`
(`DISATLAS_POWERLAW_BOOT`, `DISATLAS_KSEL_RUN`). `test_data_index.py` checks that every tracked data file is read by shipped
code and that `data/DATA_INDEX.tsv` is current (`python scripts/maintainer/data_index.py`).
