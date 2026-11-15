# Errata: manuscript numbers that the released outputs do not reproduce

Every number in the main text is listed in `tests/main_text_numbers.tsv`. `tests/test_main_parity.py` checks each
`computed` row against the released model outputs. A printed value that comes from an earlier, superseded readout
has kind `erratum`. `test_main_text_erratum` checks two things for each one:

- the corrected value (`key`) is computed from the released outputs and differs from the printed value;
- where the stale source is shipped (`stale=<key>` in the note), the stale source reproduces the printed value.

Once the text is corrected, change the row back to `computed`. The test fails if a printed value already matches
the corrected value, so a fixed row cannot stay marked as an erratum by mistake.

Every figure panel and every table is generated from the released five-seed gene-anchored cosine. Only Table S4
prints values that the released outputs do not reproduce (E4).

## Decimal tolerance

Decimal values (AUROC, percentile, share, ratio, difference) whose printed value differs from the regenerated value by
at most `DECIMAL_TOL = 0.003` (absolute, in the printed unit; `tests/tolerance.py`) are `computed` rows with op
`dtol`. They are not errata. The note of such a row starts with `tol_note:` and gives both values. The tolerance is
not used for counts, p-values, whole-percent figures or odds ratios.

The single constant moves five printed values, all checked through `dtol_ok`:

| Where | Printed | Regenerated |
|---|---|---|
| Main text, N068 (`rc_mean_pct_common_anchor`) | 0.934 | 0.93467 |
| SI registry, SI01018 (S9, `ra4357\|DisAtlas\|min`) | 0.967 | 0.96752 |
| S4.5, outside-Orphanet AUROC, all pairs (`tests/test_sibench_parity.py`) | 0.876 | 0.878 |
| S4.5, outside-Orphanet AUROC, test split | 0.909 | 0.911 |
| S4.5, HPO Lin-BMA difference, upper 95% bound | 0.122 | 0.123 |

The shuffled-label AMI ("at most 0.003"; regenerated maximum 0.002 over 50 shuffles) is an upper bound that the
regenerated value satisfies; it is checked as `le 0.003` and is neither an erratum nor a tolerance case.

## Resolved in the final text (main text)

The final revised text prints the released values for the earlier entries E1-E3, and the rows are `computed` again:

- E1, Discussion event-time odds ratio: the sentence that printed 1.93 (95% CI 1.90-1.97) is no longer in the text.
  The text gives the released 2.20 (2.15-2.25).
- E2, illustrative cosines: Crohn disease - ulcerative colitis 0.907 (N052); diabetes insipidus - type 2 diabetes
  0.181 (N053).
- E3, control-pair scores: LAL deficiency (Wolman) - CESD 0.735 at the 99.5th percentile (M122, M123); NF1 - NF2
  0.354 (M128). The superseded decoder scores are no longer quoted or shipped; `data/derived/control_pairs/control_pairs.tsv`
  lists the pairs only.


## E5. Methods: AUROC on the development pairs

| Registry | Printed | Released value |
|---|---|---|
| N281 | AUROC 0.804 on the 1,103 development pairs | 0.807 |

The released value is the AUROC of the five-seed mean cosine (`auroc_dev`), the same score definition as every other
AUROC in the text. The mean of the five per-seed AUROCs is 0.805. The gap to the printed value (0.0033) exceeds
`DECIMAL_TOL`. No shipped source reproduces 0.804.

## E4. Supplementary Information: values that the released outputs do not reproduce

Checked against the final revised SI (`disatlas_representation_learning_story.md`, sha256 recorded in `tests/README.md`).
Each printed value is an `erratum` row in `tests/si_text_numbers.tsv` whose `key` computes the released value;
`SI_BENCH_MISMATCH` in `tests/test_sibench_parity.py` is empty, because every entry of the earlier register is either
printed with its released value in the final text or no longer in the text.

| Where | Printed | Released value | Registry rows |
|---|---|---|---|
| S4.5 and Table S4, HPO-only classifier: DisAtlas gain over the four HPO measures | +0.06 to +0.08 | +0.059 to +0.086 (upper bound rounds to 0.09) | SI06794, SI06822 |
| S4.5 and Table S4: HPO measures whose 95% CI excludes zero | 3 of 4 | 4 of 4 | SI06824 |
| S4.5, HPO Lin-BMA difference, lower 95% bound | -0.001 | +0.002 (0.0017; the interval excludes zero) | SI06796 |
| S4.5 and Table S4, least-concordant pairs: DisAtlas in-graph AUROC | 0.68-0.76 | 0.686-0.766 (rounds to 0.69-0.77) | SI06808, SI06809, SI06835, SI06836 |
| S9, seed SD against disease-cluster bootstrap SD | 7.8 times | 7.95 (7.8 is the ratio of the rounded SDs 0.0180/0.0023) | SI00944 |

The S4.5 and Table S4 values were printed from the retired decoder score with an earlier Poisson bootstrap; the
released values use the cosine score with the shared disease-cluster bootstrap (B = 2,000, seed 4707). Table S4 is the
only table with corrected values (`release/coverage_table.tsv`). The S13 top-10% precision paragraph printed four DisAtlas
means from an earlier scoring run of the event-time models (0.57%, 0.37%, 9.91%, 5.50%); the author has updated them to
the released run (0.53%, 0.35%, 10.21%, 5.79%; applied by author, 2026-11-15), so registry rows SI05372, SI05375, SI05376
and SI05379 are `computed` rows again (`release/text_edits_for_author.tsv`).

**Supplementary Table S10 / SI cluster themes.** The two-sided Mann-Whitney p for gene-share fold (anatomy-defined versus
process-defined clusters) is printed as 0.032 (medians 7.8 and 17.2). It was computed on the authors' disease-gene table
(an input of the broad channel), which contains 610 pairs that cannot be redistributed. On the released table the same
test gives p = 0.052 (medians 6.1 and 13.1); the direction of the difference is unchanged, and the released-table p is
close to the 0.05 level (exact and permutation p-values differ in the third decimal). The HPO (printed
0.005) and GO-process (printed 0.001) contrasts become 0.004 and 0.002. No theme assignment of any of the 55 clusters
changes. The printed values reproduce from the authors' tables, so they stay bound to the shipped Table S10 values and
are not `erratum` rows; [data/README.md](../data/README.md) ("Table S10 (cluster themes) and the disease-gene table")
gives the input differences and the GO rank-1 match (42 of 55).

The earlier E4 entries (S5.4 ablation arms, S7 gene-anchoring AUROC, S9 placement confidence and degree controls, S9
k-medoids bootstrap, S10 node-disjoint rare pairs, negative-set realism, and the S13 event-time tables) are resolved:
the final text prints the released values, or the passage is no longer in the text. The event-time tables
(Supplementary Tables S13-S16) are regenerated row for row.

## Numbers awaiting a shipped source

The main-text registry has no `pending` rows. The SI registry has 194, by section: S2 6, S4 14, S5 2, S6 21, S7 40,
S9 27, S10 2, S11 71, S12 4, S20 7. Each `pending` row names the missing source in its note. The largest groups are
the S11 GO aggregation-rule alternatives and umbrella-term placement counts (five retrained networks and their
placement tables are not shipped), and in S7 the gene-anchoring paragraph (final model against a retrain that omits the
rule, not shipped) and the gene-anchored mechanism-class and complement-set control readouts. Pending rows are
skipped by the parity tests; `test_no_stale_pending` fails as soon as one of these values appears in a shipped value
table, so a row cannot stay pending once its source is added.

## Corrected values overlaid on the printed text

`release/manuscript_overlay.tsv` lists every printed manuscript or Supplementary value that the correction notes 01-04
replace (new note cycle, numbered from 01): label, location, kind, old (printed) value, new value and note. These values
replace the printed values and are the reproduction target for DisAtlas. The table parity check compares the regenerated
tables with the printed tables after the overlay is applied; `tests/test_manuscript_overlay.py` checks that every new
table row is regenerated and that no old row remains.
