# IxIDN-Matched benchmark

- `ixidn_matched_4364.tsv` — the IxIDN-Matched benchmark: 4,364 pairs (2,182 IxIDN+ / 2,182 IxIDN-)
  over 938 diseases, with the tune/test pool of each pair. `code/01_build_ixidn/04_match_ixidn.py`
  rebuilds it as `data/01_rebuilt/ixidn_matched.tsv`; `tests/test_ixidn_parity.py` checks the two
  are identical.
