"""Appendix, input data sources -- Table tbl:SI_input_dis_symptom_sources.

Per-source breakdown of the tier-1 disease-symptom product (after capping), from
data/00_raw/associations/disease_symptom_tier1_sources.parquet; the MONDO-propagation row is
the propagated evidence tier (evidence_tier == 9) of data/00_raw/associations/disease_symptom.parquet.
Percentages are of the 13,080-disease candidate universe.
"""
import sys
from pathlib import Path

import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _sources_common import ASSOC, candidate_universe, per_source  # noqa: E402
from tablib import thousands, pct1, write_and_check  # noqa: E402

U = candidate_universe()
assert len(U) == 13080, len(U)
rows, comb = per_source("disease_symptom", "hpo_id")
out = [[r["label"], thousands(r["assoc"]), thousands(r["diseases"]), thousands(r["items"]),
        pct1(r["diseases"] / len(U))] for r in rows]
out.append(["**Combined (tier-1, after capping)**", f"**{thousands(comb['assoc'])}**",
            f"**{thousands(comb['diseases'])}**", f"**{thousands(comb['items'])}**",
            f"**{pct1(comb['diseases'] / len(U))}**"])
full = pd.read_parquet(ASSOC / "disease_symptom.parquet", columns=["disease_id", "evidence_tier"])
prop = full[full.evidence_tier == 9]
out.append(["MONDO-propagation enrichment (secondary augmentation, not an independent source)",
            thousands(len(prop)), thousands(prop.disease_id.nunique()), "—",
            pct1(prop.disease_id.nunique() / len(U))])
df = pd.DataFrame(out, columns=["Source", "Associations", "Diseases", "Symptoms", "% of 13,080 candidates"])

EXPECTED = [
    ["BIKG (HPOA-derived)", "147,952", "7,031", "9,622", "53.8%"],
    ["OMIM (phenotype HPOA)", "123,388", "8,099", "8,800", "61.9%"],
    ["**Combined (tier-1, after capping)**", "**167,798**", "**8,991**", "**10,199**", "**68.7%**"],
    ["MONDO-propagation enrichment (secondary augmentation, not an independent source)",
     "39,894", "986", "—", "7.5%"],
]

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_SI_input_dis_symptom_sources", EXPECTED,
                    "Disease-symptom association sources feeding the DisAtlas pipeline and the MONDO "
                    "propagation layer. Percentages are of the pipeline's 13,080-candidate universe.")
