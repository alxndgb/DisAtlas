# figD: parity notes

Figures: release/fig_parity_figD.tsv. Tables: release/table_parity_figD.tsv. Prose numbers: release/numbers_figD.tsv
(20 of 20 text-checkable values occur verbatim in the final manuscript; three are figure-panel values, checked through the
byte-identical figure; one is the decoder-era S13 subset row, listed as a residual).

## Test suite (pytest tests/)

Branch base: 26 failed, 17 errors. This branch: 30 failed, 17 errors. The 43 shared failures (test_main_parity.py,
test_si_event_time_or.py) are unchanged by this branch.

The 4 failures that are new to this branch come from fixtures in tests/parity/sibio_event_time.json and sibio_fda.json:
- SI_fda_repurposing, SI_fda_prioritisation, SI_fda_prioritisation_text: the fixture strings come from an earlier manuscript
  draft (for example 91.9% and -6.2 percentile points). Neither occurs in the final manuscript. The final text gives 90.7%
  (87.2-93.8%), 63.9% and -7.2 percentile points, and these equal the regenerated gene-anchored outputs.
- SI_et_coverage_text: the fixture holds the decoder-era Supplementary Table S10 body, which HEAD carried as a stale render.
  The regenerated table follows the shipped gene-anchored deltas.tsv (see the S10/S13 rows of table_parity_figD.tsv).

These fixtures belong to the shared SI-bio parity suite, so this branch leaves them unchanged and the integration branch can
decide on them. test_sibench_parity.py has been updated to the gene-anchored values in the final text and passes in full
(24/24).
