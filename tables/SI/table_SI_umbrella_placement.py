"""Appendix, "Placement of umbrella disease terms" -- text values, and the released per-disease flag.

Inputs:
  tables/inputs/umbrella_placement.tsv     the 444 umbrella terms (>= 5 lexical subtypes) of the network:
                                           subtype_average_other_cluster = the subtype average differs from the
                                           disease's position by more than projection noise and lies in another
                                           cluster (93); alternative_position_flag = additionally no own symptom
                                           data (69), with the subtype-average cluster, shift and noise
  tables/inputs/cluster_labels_by_K.tsv    K = 50 cluster of every network disease
Recomputed: 444, 93, 69, 24 = 93 - 69, 75% = (444 - 112) / 444, flagged shift > noise, flagged clusters differ.
Stated from the screen and retraining runs (the MONDO hierarchy and retrained models are not shipped):
  112 shifted terms, 10 of 93 explained by the inherited-symptom cap, 613 capped diseases, 111 restored at a cap
  of 30 descendants, -51% case-study score (82 times the noise floor, 4 of 5 seeds), 537 of 2,499 grouping-class
  diseases without own symptom annotation.
"""
import sys
from pathlib import Path

import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import INPUTS, to_markdown  # noqa: E402

U = pd.read_csv(INPUTS / "umbrella_placement.tsv", sep="\t")
K50 = pd.read_csv(INPUTS / "cluster_labels_by_K.tsv", sep="\t").set_index("disease_id").K50
STATED = {"shifted": 112, "cap_explained": 10, "capped": 613, "cap30_restored": 111}

assert len(U) == U.disease_id.nunique() == 444 and U.disease_id.isin(K50.index).all()
assert (U.n_lexical_subtypes >= 5).all()
other, flag = U[U.subtype_average_other_cluster == 1], U[U.alternative_position_flag == 1]
assert set(flag.disease_id) <= set(other.disease_id)
assert (flag.map_shift > flag.projection_noise).all()
assert (flag.subtype_average_cluster.astype(int).values != flag.disease_id.map(K50).values).all()
assert U.subtype_average_cluster.notna().sum() == len(flag)
noise = (len(U) - STATED["shifted"]) / len(U)

text = pd.DataFrame([
    ["umbrella terms with >= 5 subtypes", f"{len(U)}", "444"],
    ["shifted beyond projection noise", f"{STATED['shifted']} ({100 * STATED['shifted'] / len(U):.1f}%)", "112 (25.2%)"],
    ["shifted and in a different cluster", f"{len(other)}", "93"],
    ["of these, explained by the inherited-symptom cap", f"{STATED['cap_explained']} of {len(other)}", "10 of the 93"],
    ["diseases under the cap / restored at 30 descendants", f"{STATED['capped']} / {STATED['cap30_restored']}", "613 / 111"],
    ["case-study score change at cap 30", "-51% (82 x noise floor, 4 of 5 seeds)", "51% lower (82 times the noise floor, 4 of 5 seeds)"],
    ["umbrella terms with no real signal", f"{100 * noise:.0f}%", "75%"],
    ["mismatched diseases with own symptom data (would be overwritten)", f"{len(other) - len(flag)} of {len(other)}", "24 of the 93"],
    ["flagged with a subtype-based alternative position", f"{len(flag)}", "69"],
    ["grouping-class diseases without own symptom annotation", "537 of 2,499", "537 of the 2,499"],
], columns=["quantity", "rebuilt", "manuscript"])
EXPECTED_REBUILT = ["444", "112 (25.2%)", "93", "10 of 93", "613 / 111", "-51% (82 x noise floor, 4 of 5 seeds)",
                    "75%", "24 of 93", "69", "537 of 2,499"]
assert text.rebuilt.tolist() == EXPECTED_REBUILT, text.rebuilt.tolist()

if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    text.to_csv(here / "tbl_SI_umbrella_placement_text_values.csv", index=False)
    (here / "tbl_SI_umbrella_placement_text_values.md").write_text(to_markdown(text))
    print(text.to_string(index=False))
