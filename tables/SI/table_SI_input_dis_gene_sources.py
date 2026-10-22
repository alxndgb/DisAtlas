"""Appendix, input data sources -- Table tbl:SI_input_dis_gene_sources.

Per-source breakdown of the tier-1 disease-gene product (after capping), from
data/00_raw/associations/disease_gene_tier1_sources.parquet. Percentages are of the
13,080-disease candidate universe (see _sources_common.candidate_universe).

Note: the column printed as "Diseases exclusively from this source" counts tier-1
disease-gene ASSOCIATIONS supported by that source alone (the quantity reported in the
manuscript); see tables/README.md.
"""
import sys
from pathlib import Path

import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from _sources_common import candidate_universe, per_source  # noqa: E402
from tablib import thousands, pct1, write_and_check  # noqa: E402

U = candidate_universe()
assert len(U) == 13080, len(U)
rows, comb = per_source("disease_gene", "gene_id")
out = [[r["label"], thousands(r["assoc"]), thousands(r["diseases"]), thousands(r["items"]),
        thousands(r["solo_assoc"]), pct1(r["diseases"] / len(U))] for r in rows]
out.append(["**Combined (tier-1, after capping)**", f"**{thousands(comb['assoc'])}**",
            f"**{thousands(comb['diseases'])}**", f"**{thousands(comb['items'])}**", "—",
            f"**{pct1(comb['diseases'] / len(U))}**"])
df = pd.DataFrame(out, columns=["Source", "Associations", "Diseases", "Genes",
                                "Diseases exclusively from this source", "% of 13,080 candidates"])

EXPECTED = [
    ["OpenTargets", "8,009", "5,421", "4,063", "1,693", "41.4%"],
    ["GenCC", "6,903", "6,033", "4,157", "1,232", "46.1%"],
    ["OMIM (direct)", "6,019", "5,818", "4,328", "588", "44.5%"],
    ["Orphanet", "8,150", "4,090", "4,421", "4,183", "31.3%"],
    ["BIKG (NLP-mined scientific literature)", "3,413", "783", "995", "2,818", "6.0%"],
    ["**Combined (tier-1, after capping)**", "**18,287**", "**9,426**", "**5,656**", "—", "**72.1%**"],
]

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_SI_input_dis_gene_sources", EXPECTED,
                    "Disease-gene association sources feeding the DisAtlas pipeline. "
                    "Percentages are of the pipeline's 13,080-candidate universe.")
