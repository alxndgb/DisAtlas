# Figure-specific inputs (figures/main/data)

Inputs used only by the main figures and not shipped elsewhere in the repository.

- `fig5_display_map.tsv.gz`: `disease_id`, `cluster`, `x`, `y` for the 13,019 network diseases.
  `cluster` = K = 50 UPNA-DDAL cluster; `x`, `y` = weighted-MDS display coordinates (Figure 5a/b).
- `fig5_control_pairs.tsv`: the five validated control pairs in Figure 5a (`role` pos/neg) and their
  UPNA-DDAL probability `p` from the full probability matrix (ticks in Figure 5d).
- `ixidn_matched_go_resnik.tsv.gz`: `blind_id`, `go_resnik_genebma` (GO Resnik gene-BMA similarity per
  IxIDN-Matched pair; blank = not scorable). Baseline in Figure 5c.
