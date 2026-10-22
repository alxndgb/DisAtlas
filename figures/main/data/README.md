# Figure-specific inputs (figures/main/data)

Inputs used only by the main figures and not shipped elsewhere in the repository.

- `fig3_display_map.tsv.gz`: `disease_id`, `cluster`, `x`, `y` for the 13,019 network diseases.
  `cluster` = K = 50 UPNA-DDAL cluster; `x`, `y` = weighted-MDS display coordinates (Figure 3a/b).
- `fig3_control_pairs.tsv`: the five validated control pairs in Figure 3a (`role` pos/neg) and their
  UPNA-DDAL probability `p` from the full probability matrix (ticks in Figure 3d).
- `ixidn_matched_go_resnik.tsv.gz`: `blind_id`, `go_resnik_genebma` (GO Resnik gene-BMA similarity per
  IxIDN-Matched pair; blank = not scorable). Baseline in Figure 3c.
