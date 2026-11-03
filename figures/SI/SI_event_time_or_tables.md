| Decision year | Evidence up to | MS rank | MS top % | NMO rank | NMO top % |
|---:|---:|---:|---:|---:|---:|
| 1996 | 1995 | 32 of 644 | 5.0 | not in map |  |
| 1997 | 1996 | 39 of 741 | 5.3 | not in map |  |
| 1999 | 1998 | 47 of 923 | 5.1 | not in map |  |
| 2000 | 1999 | 112 of 1,041 | 10.8 | not in map |  |
| 2001 | 2000 | 126 of 1,155 | 10.9 | not in map |  |
| 2002 | 2001 | 200 of 1,260 | 15.9 | not in map |  |
| 2003 | 2002 | 289 of 1,367 | 21.1 | not in map |  |
| 2004 | 2003 | 259 of 1,483 | 17.5 | not in map |  |
| 2005 | 2004 | 294 of 1,576 | 18.7 | not in map |  |
| 2006 | 2005 | 101 of 1,648 | 6.1 | not in map |  |
| 2007 | 2006 | 108 of 1,728 | 6.2 | 45 of 1,728 | 2.6 |
| 2008 | 2007 | 312 of 1,809 | 17.2 | 87 of 1,809 | 4.8 |
| 2009 | 2008 | 415 of 1,874 | 22.1 | 97 of 1,874 | 5.2 |

: Eculizumab, DisAtlas rank of MS and NMO among diseases held out of benchmark tuning, by decision year.
Each decision year uses the map built from evidence up to the year before.
Rank 1 is the closest disease to PNH, and top % is the rank as a share of candidates.
Not in map means that NMO has no node in that snapshot.
There is no row for 1998. {#tbl:SI_nmo_ms_by_slice}

| Decision year | New trials | Drugs | DisAtlas, odds ratio per 1 SD (95% CI) | Weighted Adamic–Adar, odds ratio per 1 SD (95% CI) |
|---|---|---|---|---|
| 2000 | 675 | 352 | 1.87 (1.71–2.06) | 1.12 (1.08–1.17) |
| 2001 | 752 | 389 | 1.73 (1.59–1.89) | 1.18 (1.14–1.22) |
| 2002 | 900 | 450 | 1.65 (1.50–1.82) | 1.12 (1.08–1.16) |
| 2003 | 914 | 516 | 1.87 (1.73–2.03) | 1.15 (1.11–1.20) |
| 2004 | 1,117 | 585 | 1.95 (1.79–2.13) | 1.14 (1.10–1.18) |
| 2005 | 1,283 | 616 | 1.89 (1.73–2.06) | 1.09 (1.06–1.12) |
| 2006 | 1,277 | 654 | 1.72 (1.58–1.86) | 1.10 (1.07–1.13) |
| 2007 | 1,308 | 723 | 1.99 (1.84–2.15) | 1.09 (1.06–1.11) |
| 2008 | 1,316 | 731 | 2.23 (2.05–2.42) | 1.13 (1.10–1.16) |
| 2009 | 1,558 | 832 | 1.97 (1.83–2.12) | 1.11 (1.09–1.14) |
| 2010 | 1,543 | 829 | 2.09 (1.94–2.24) | 1.08 (1.06–1.11) |
| 2011 | 1,817 | 865 | 1.79 (1.65–1.94) | 1.11 (1.09–1.13) |
| 2012 | 1,500 | 847 | 2.20 (2.03–2.38) | 1.13 (1.10–1.16) |
| 2013 | 1,410 | 832 | 1.94 (1.78–2.12) | 1.11 (1.07–1.14) |
| 2014 | 1,656 | 926 | 1.89 (1.73–2.06) | 1.12 (1.09–1.16) |
| 2015 | 1,937 | 978 | 2.23 (2.06–2.40) | 1.14 (1.12–1.17) |
| 2016 | 2,173 | 1,056 | 2.06 (1.93–2.20) | 1.09 (1.07–1.11) |
| 2017 | 1,977 | 1,027 | 1.97 (1.85–2.11) | 1.10 (1.08–1.12) |
| 2018 | 2,131 | 1,138 | 2.21 (2.04–2.38) | 1.11 (1.09–1.13) |
| 2019 | 2,193 | 1,166 | 2.16 (2.01–2.33) | 1.12 (1.10–1.14) |
| 2020 | 2,708 | 1,348 | 1.62 (1.53–1.72) | 1.05 (1.03–1.06) |
| 2021 | 2,447 | 1,374 | 1.98 (1.86–2.11) | 1.09 (1.07–1.10) |
| 2000–2021 pooled | 34,592 | 6,469 | 1.93 (1.90–1.97) | 1.10 (1.10–1.11) |

: Event-time holdout, odds that a disease not yet trialed with a drug is the one it is trialed in next, per 1 SD of score, each adjusted for popularity, by decision year and pooled, among diseases held out of benchmark tuning.
DisAtlas values are the dark blue series of [Figure @fig:event_time]b; weighted Adamic–Adar is computed on the same yearly MechDist graph without learning.
Conditional logistic model with one stratum per drug and decision year (Efron ties), popularity as a natural spline plus a zero indicator; 95% intervals from robust standard errors clustered by drug.
Decision years share drugs, so the yearly estimates are not independent. {#tbl:SI_event_time_or_yearly}

| Method | Comparison set | Candidate diseases | Probability (95% CI) | DisAtlas minus method (95% CI) |
|---|---|---|---|---|
| DisAtlas | Popularity- and connectivity-matched | Held out of benchmark tuning | 0.682 (0.677–0.686) | — |
|  | Popularity- and connectivity-matched | All diseases | 0.685 (0.681–0.689) | — |
|  | Popularity-matched | Held out of benchmark tuning | 0.680 (0.676–0.685) | — |
|  | Popularity-matched | All diseases | 0.690 (0.686–0.694) | — |
| Weighted Adamic–Adar (same graph, no learning) | Popularity- and connectivity-matched | Held out of benchmark tuning | 0.641 (0.636–0.645) | 0.041 (0.036–0.045) |
|  | Popularity- and connectivity-matched | All diseases | 0.665 (0.661–0.669) | 0.021 (0.017–0.024) |
|  | Popularity-matched | Held out of benchmark tuning | 0.627 (0.623–0.631) | 0.053 (0.049–0.058) |
|  | Popularity-matched | All diseases | 0.643 (0.639–0.647) | 0.047 (0.044–0.050) |
| Disease popularity | Popularity- and connectivity-matched | Held out of benchmark tuning | 0.619 (0.614–0.623) | 0.063 (0.056–0.069) |
|  | Popularity- and connectivity-matched | All diseases | 0.601 (0.598–0.605) | 0.084 (0.078–0.089) |
|  | Popularity-matched | Held out of benchmark tuning | 0.621 (0.617–0.625) | 0.059 (0.053–0.066) |
|  | Popularity-matched | All diseases | 0.609 (0.605–0.613) | 0.081 (0.076–0.086) |
| Shared genes | Popularity- and connectivity-matched | Held out of benchmark tuning | 0.647 (0.643–0.652) | 0.034 (0.030–0.038) |
|  | Popularity- and connectivity-matched | All diseases | 0.664 (0.660–0.667) | 0.021 (0.018–0.025) |
|  | Popularity-matched | Held out of benchmark tuning | 0.650 (0.646–0.655) | 0.030 (0.026–0.034) |
|  | Popularity-matched | All diseases | 0.670 (0.667–0.674) | 0.020 (0.017–0.023) |
| GO biological-process similarity | Popularity- and connectivity-matched | Held out of benchmark tuning | 0.619 (0.615–0.623) | 0.062 (0.058–0.067) |
|  | Popularity- and connectivity-matched | All diseases | 0.621 (0.617–0.624) | 0.065 (0.061–0.068) |
|  | Popularity-matched | Held out of benchmark tuning | 0.619 (0.615–0.623) | 0.061 (0.057–0.065) |
|  | Popularity-matched | All diseases | 0.628 (0.624–0.631) | 0.062 (0.059–0.066) |
| HPO phenotype similarity | Popularity- and connectivity-matched | Held out of benchmark tuning | 0.515 (0.512–0.519) | 0.166 (0.161–0.171) |
|  | Popularity- and connectivity-matched | All diseases | 0.523 (0.520–0.526) | 0.162 (0.158–0.167) |
|  | Popularity-matched | Held out of benchmark tuning | 0.515 (0.512–0.519) | 0.165 (0.160–0.170) |
|  | Popularity-matched | All diseases | 0.523 (0.520–0.526) | 0.167 (0.163–0.171) |

: Event-time holdout, all decision years pooled, on the per-drug probability scale: probability that each method scores the next trialed disease above a matched disease not yet trialed with the drug (random ranking, 0.5), averaged over drug–decision-year pairs, and the paired difference from DisAtlas.
Brackets give 95% intervals from a bootstrap over drugs (2,000 replicates, the replicates of the published analysis).
The pre-specified test is the median-percentile comparison with weighted Adamic–Adar ([Supplementary Table @tbl:SI_event_time_deltas]); this table re-expresses it and the other reference scores.
Shared-gene, GO biological-process and HPO phenotype similarity use current annotation releases for every decision year, so unlike the other scores they are not limited to earlier evidence. {#tbl:SI_event_time_auc}

| Decision year | Drugs | Popularity- and connectivity-matched: mean (95% CI) | Popularity- and connectivity-matched: median (interquartile range) | Popularity-matched: mean (95% CI) | Popularity-matched: median (interquartile range) | No matching: mean (95% CI) | No matching: median (interquartile range) |
|---|---|---|---|---|---|---|---|
| 1996 | 117 | 0.607 (0.557–0.657) | 0.637 [0.392–0.833] | 0.609 (0.559–0.658) | 0.637 [0.400–0.833] | 0.800 (0.758–0.836) | 0.871 [0.720–0.941] |
| 1997 | 167 | 0.624 (0.583–0.661) | 0.625 [0.437–0.854] | 0.624 (0.583–0.662) | 0.625 [0.437–0.854] | 0.805 (0.768–0.836) | 0.879 [0.777–0.945] |
| 1999 | 295 | 0.624 (0.594–0.652) | 0.650 [0.447–0.833] | 0.625 (0.595–0.654) | 0.651 [0.441–0.833] | 0.825 (0.801–0.847) | 0.896 [0.792–0.949] |
| 2000 | 352 | 0.704 (0.681–0.727) | 0.734 [0.550–0.889] | 0.704 (0.681–0.727) | 0.732 [0.556–0.889] | 0.862 (0.844–0.879) | 0.924 [0.836–0.968] |
| 2001 | 389 | 0.665 (0.639–0.692) | 0.717 [0.481–0.889] | 0.664 (0.639–0.691) | 0.717 [0.481–0.889] | 0.813 (0.791–0.834) | 0.906 [0.750–0.971] |
| 2002 | 450 | 0.667 (0.644–0.689) | 0.682 [0.500–0.870] | 0.668 (0.645–0.690) | 0.684 [0.500–0.870] | 0.833 (0.814–0.851) | 0.906 [0.786–0.968] |
| 2003 | 516 | 0.673 (0.649–0.695) | 0.726 [0.506–0.894] | 0.673 (0.649–0.695) | 0.720 [0.500–0.894] | 0.852 (0.836–0.868) | 0.921 [0.812–0.977] |
| 2004 | 585 | 0.678 (0.657–0.700) | 0.733 [0.500–0.894] | 0.678 (0.656–0.699) | 0.733 [0.500–0.889] | 0.862 (0.847–0.876) | 0.933 [0.830–0.975] |
| 2005 | 616 | 0.684 (0.662–0.704) | 0.731 [0.509–0.915] | 0.684 (0.662–0.704) | 0.735 [0.509–0.913] | 0.861 (0.846–0.876) | 0.931 [0.822–0.981] |
| 2006 | 654 | 0.655 (0.633–0.675) | 0.696 [0.463–0.901] | 0.650 (0.628–0.671) | 0.692 [0.455–0.900] | 0.843 (0.827–0.859) | 0.917 [0.798–0.979] |
| 2007 | 723 | 0.681 (0.664–0.701) | 0.748 [0.500–0.900] | 0.680 (0.663–0.700) | 0.741 [0.500–0.899] | 0.852 (0.839–0.866) | 0.932 [0.810–0.976] |
| 2008 | 731 | 0.685 (0.667–0.704) | 0.739 [0.536–0.902] | 0.682 (0.663–0.700) | 0.742 [0.527–0.900] | 0.872 (0.860–0.885) | 0.942 [0.853–0.974] |
| 2009 | 832 | 0.666 (0.648–0.684) | 0.715 [0.500–0.897] | 0.661 (0.643–0.679) | 0.703 [0.492–0.881] | 0.860 (0.848–0.872) | 0.923 [0.814–0.973] |
| 2010 | 829 | 0.681 (0.664–0.698) | 0.740 [0.536–0.888] | 0.679 (0.663–0.696) | 0.744 [0.526–0.890] | 0.871 (0.860–0.881) | 0.925 [0.845–0.970] |
| 2011 | 865 | 0.653 (0.635–0.670) | 0.694 [0.500–0.864] | 0.651 (0.633–0.669) | 0.696 [0.486–0.872] | 0.874 (0.863–0.884) | 0.926 [0.842–0.969] |
| 2012 | 847 | 0.679 (0.663–0.696) | 0.724 [0.521–0.893] | 0.691 (0.675–0.707) | 0.744 [0.529–0.902] | 0.892 (0.883–0.901) | 0.941 [0.867–0.978] |
| 2013 | 832 | 0.662 (0.643–0.680) | 0.716 [0.500–0.875] | 0.667 (0.648–0.684) | 0.708 [0.500–0.884] | 0.887 (0.877–0.897) | 0.933 [0.859–0.976] |
| 2014 | 926 | 0.669 (0.653–0.686) | 0.714 [0.493–0.886] | 0.670 (0.655–0.687) | 0.714 [0.489–0.886] | 0.887 (0.878–0.896) | 0.939 [0.867–0.979] |
| 2015 | 978 | 0.709 (0.693–0.724) | 0.760 [0.563–0.923] | 0.713 (0.697–0.727) | 0.765 [0.561–0.920] | 0.899 (0.889–0.907) | 0.944 [0.872–0.983] |
| 2016 | 1,056 | 0.712 (0.697–0.726) | 0.789 [0.563–0.912] | 0.698 (0.683–0.713) | 0.764 [0.543–0.909] | 0.908 (0.901–0.915) | 0.949 [0.889–0.983] |
| 2017 | 1,027 | 0.700 (0.686–0.715) | 0.767 [0.561–0.909] | 0.689 (0.674–0.704) | 0.739 [0.537–0.902] | 0.902 (0.894–0.910) | 0.949 [0.884–0.980] |
| 2018 | 1,138 | 0.712 (0.697–0.726) | 0.778 [0.564–0.909] | 0.708 (0.695–0.722) | 0.759 [0.562–0.906] | 0.916 (0.910–0.923) | 0.957 [0.901–0.984] |
| 2019 | 1,166 | 0.710 (0.696–0.725) | 0.786 [0.570–0.917] | 0.714 (0.701–0.728) | 0.761 [0.566–0.914] | 0.918 (0.912–0.925) | 0.957 [0.901–0.986] |
| 2020 | 1,348 | 0.661 (0.646–0.675) | 0.719 [0.500–0.885] | 0.648 (0.634–0.663) | 0.693 [0.457–0.879] | 0.892 (0.885–0.900) | 0.941 [0.862–0.982] |
| 2021 | 1,374 | 0.689 (0.675–0.703) | 0.758 [0.544–0.895] | 0.700 (0.688–0.713) | 0.756 [0.558–0.898] | 0.918 (0.912–0.925) | 0.961 [0.910–0.986] |
| All years | 18,813 pairs | 0.682 (0.677–0.686) | — | 0.680 (0.676–0.685) | — | 0.882 (0.880–0.885) | — |

: Event-time holdout, DisAtlas, by decision year and pooled: per-drug probability that the next trialed disease is scored above an alternative disease not yet trialed with the drug, in each comparison set, in diseases held out of benchmark tuning, in the pre-registered matched analysis.
Mean over drugs with its 95% bootstrap interval over drugs (2,000 replicates), and the median and interquartile range of the per-drug values; the last row averages over all drug–decision-year pairs.
Drugs are counted in the popularity- and connectivity-matched set; a drug whose only new trials of the year had no matched alternative is not counted in that set.
There is no row for 1998. {#tbl:SI_event_time_auc_yearly}

| Decision year | Evidence up to | Drugs | New trials | New trials per drug | Comparison set | Median candidate diseases | Trials with fewer than 10 candidates (%) | Drugs with a new trial in the top 10% (%, 95% CI) | Random (%) | Probability (95% CI) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1996 | 1995 | 117 | 240 | 2.05 | Popularity-matched | 9 | 52.1 | 11.1 (5.1–17.1) | 7.9 | 0.609 (0.561–0.656) |
| | | | | | Popularity- and connectivity-matched | 9 | 52.1 | 11.1 (6.0–17.1) | 8.0 | 0.607 (0.557–0.655) |
| | | | | | No matching | 640 | 0.0 | 61.5 (53.0–70.1) | 18.0 | 0.800 (0.761–0.835) |
| 1997 | 1996 | 168 | 315 | 1.88 | Popularity-matched | 12 | 19.4 | 26.8 (20.2–33.3) | 11.6 | 0.624 (0.584–0.664) |
| | | | | | Popularity- and connectivity-matched | 12 | 19.4 | 26.8 (20.2–33.3) | 11.8 | 0.624 (0.583–0.664) |
| | | | | | No matching | 737 | 0.0 | 66.7 (59.5–73.2) | 17.0 | 0.805 (0.770–0.838) |
| 1999 | 1998 | 295 | 646 | 2.19 | Popularity-matched | 15 | 17.8 | 19.7 (15.3–24.1) | 11.7 | 0.625 (0.595–0.655) |
| | | | | | Popularity- and connectivity-matched | 15 | 17.8 | 20.0 (15.6–24.7) | 12.1 | 0.624 (0.593–0.653) |
| | | | | | No matching | 918 | 0.0 | 71.2 (65.8–76.6) | 19.2 | 0.825 (0.801–0.847) |
| 2000 | 1999 | 352 | 675 | 1.92 | Popularity-matched | 17 | 4.3 | 23.9 (19.6–28.7) | 12.4 | 0.704 (0.681–0.728) |
| | | | | | Popularity- and connectivity-matched | 17 | 4.3 | 23.9 (19.6–28.4) | 12.6 | 0.704 (0.681–0.727) |
| | | | | | No matching | 1036 | 0.0 | 73.3 (68.5–77.8) | 17.2 | 0.862 (0.844–0.880) |

: Event-time holdout, early decision years, DisAtlas, diseases held out of benchmark tuning. [Figure @fig:event_time]b plots decision years from 2000 (its first year is shown for comparison); all pooled numbers use the pre-specified decision years 1996, 1997 and 1999–2021.
Top-10% share, drugs whose best-ranked new trial of the year is in the top 10% of candidates; random, the same share expected under random ranking given each drug's number of new trials and the sizes of its matched groups; probability, the per-drug probability that the new trial is scored above a matched alternative (random, 0.5).
Brackets give 95% intervals from a bootstrap over drugs (2,000 replicates).
A trial whose matched group has fewer than 10 candidates cannot reach the top 10% (all years in [Supplementary Table @tbl:SI_event_time_small_sets]).
There is no row for 1998.
No evidence snapshot ending in 1997 was built, because snapshots before 2000 were built only at the cut-off years of the mechanisms under study and no mechanism had a 1997 cut-off, so decision year 1998 is not scored. {#tbl:SI_event_time_early_years}

| Decision year | New trials | Median matched-set size, popularity- and connectivity-matched | Median matched-set size, popularity-matched | Trials with fewer than 10 candidates (n, %) | Trials with no alternative (n) |
|---|---|---|---|---|---|
| 1996 | 240 | 9 | 9 | 125 (52.1) | 3 |
| 1997 | 315 | 12 | 12 | 61 (19.4) | 5 |
| 1999 | 646 | 15 | 15 | 115 (17.8) | 1 |
| 2000 | 675 | 17 | 17 | 29 (4.3) | 0 |
| 2001 | 752 | 19 | 19 | 49 (6.5) | 0 |
| 2002 | 900 | 21 | 22 | 51 (5.7) | 0 |
| 2003 | 914 | 24 | 24 | 33 (3.6) | 0 |
| 2004 | 1,117 | 24 | 25 | 9 (0.8) | 0 |
| 2005 | 1,283 | 24 | 24 | 32 (2.5) | 0 |
| 2006 | 1,277 | 25 | 26 | 19 (1.5) | 0 |
| 2007 | 1,308 | 23.5 | 28 | 9 (0.7) | 0 |
| 2008 | 1,316 | 25 | 31 | 7 (0.5) | 0 |
| 2009 | 1,558 | 26 | 32 | 14 (0.9) | 0 |
| 2010 | 1,543 | 26 | 34 | 9 (0.6) | 0 |
| 2011 | 1,817 | 27 | 35 | 7 (0.4) | 0 |
| 2012 | 1,500 | 26 | 37 | 18 (1.2) | 0 |
| 2013 | 1,410 | 27 | 40 | 5 (0.4) | 0 |
| 2014 | 1,656 | 28 | 43 | 3 (0.2) | 0 |
| 2015 | 1,937 | 31 | 50 | 2 (0.1) | 0 |
| 2016 | 2,173 | 31 | 53 | 1 (0.0) | 0 |
| 2017 | 1,977 | 30 | 58 | 0 (0.0) | 0 |
| 2018 | 2,131 | 28 | 63 | 0 (0.0) | 0 |
| 2019 | 2,193 | 25 | 72 | 0 (0.0) | 0 |
| 2020 | 2,708 | 27 | 77 | 0 (0.0) | 0 |
| 2021 | 2,447 | 29 | 82 | 0 (0.0) | 0 |
| 2000–2021 | 34,592 | — | — | 297 (0.9) | 0 |

: Event-time holdout, size of the matched sets by decision year, diseases held out of benchmark tuning.
The matched set of a new trial is the trialed disease and its matched alternatives not yet trialed with the drug.
A trial with fewer than 10 candidates cannot reach the top 10% of its ranking; the count is the same in both comparison sets.
Trials with no alternative are left out of all analyses.
The last row covers the decision years plotted in [Figure @fig:event_time]b.
There is no row for 1998. {#tbl:SI_event_time_small_sets}
