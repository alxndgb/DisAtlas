| Setting | Ensemble AUROC | Ensemble AUPRC | Per-seed AUROC (mean ± SD) | ΔAUROC | ΔAUPRC |
|---|---|---|---|---|---|
| DisAtlas (window 2, uniform decoder negatives, unweighted loss) | 0.782 | 0.772 | 0.776 ± 0.003 | +0.000 | +0.000 |
| Same settings, retrained (run-to-run reference) | 0.784 | 0.772 | 0.778 ± 0.007 | +0.002 | +0.001 |
| Node2vec window 1 | 0.793 | 0.778 | 0.781 ± 0.007 | +0.004 | +0.002 |
| Node2vec window 5 | 0.769 | 0.761 | 0.765 ± 0.005 | -0.012 | -0.009 |
| Degree-proportional decoder negatives | 0.784 | 0.773 | 0.778 ± 0.005 | +0.002 | +0.002 |

Sensitivity of test-split performance (1,094 pairs, 491 positive) to the two settings fixed after test-split results had been seen.
