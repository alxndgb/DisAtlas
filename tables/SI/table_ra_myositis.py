"""Supplementary Table SI_ra_myositis: rheumatoid arthritis (MONDO:0008383) against the inflammatory myopathies and
overlap disorders of the 13,021-disease index.

Mean cosine over the five seeds of the final model, range across seeds, and percentile of the mean among all
84,766,710 disease pairs; rows are ordered by mean cosine. Every disease of the group that is in the index is listed.

Input:  figures/SI/SI_cluster_rarecommon_numbers.json (key ra_case; written by figures/SI/si_cluster_rarecommon_numbers.py)
Output: table_ra_myositis.tsv and table_ra_myositis.md next to this script
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
N = json.load(open(HERE.parents[1] / "figures/SI/SI_cluster_rarecommon_numbers.json"))
RA = N["ra_case"]
# (key in ra_case, display name, MONDO identifier)
GROUP = [("dermatomyositis_mondo", "Dermatomyositis", "MONDO:0016367"), ("myositis_disease", "Myositis disease", "MONDO:0021167"),
         ("polymyositis", "Polymyositis", "MONDO:0019127"), ("antisynthetase", "Antisynthetase syndrome", "MONDO:0019344"),
         ("ibm", "Inclusion body myositis", "MONDO:0007827"), ("imnm", "Immune-mediated necrotizing myopathy", "MONDO:0016098"),
         ("cadm", "Clinically amyopathic dermatomyositis", "MONDO:1010096"), ("mctd", "Mixed connective tissue disease", "MONDO:0005854")]
rows = []
for k, name, mondo in GROUP:
    v = RA[k]
    assert v is not None, f"{k} not in the 13,021-disease index"
    rows.append((float(v["mean"]), [name, mondo, f"{float(v['mean']):.3f}", f"{float(v['lo']):.3f}–{float(v['hi']):.3f}",
                                    f"{float(v['pct']):.1f}"]))
rows = [r for _, r in sorted(rows, key=lambda x: -x[0])]
H = ["Disease", "MONDO identifier", "Mean cosine with RA", "Seed range", "Percentile among all pairs"]
pd.DataFrame(rows, columns=H).to_csv(HERE / "table_ra_myositis.tsv", sep="\t", index=False)
md = ["| " + " | ".join(H) + " |", "|---|---|---|---|---|"] + ["| " + " | ".join(r) + " |" for r in rows]
(HERE / "table_ra_myositis.md").write_text("\n".join(md) + "\n")
print(len(rows))
