"""Turn saved LLM responses into the shipped labels and check that they match.

    python code/02_llm_annotations/parse_responses.py

Reads only data/00_raw/llm_annotations/. Exits non-zero on any mismatch.
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data" / "00_raw" / "llm_annotations"
KEEP = {"STRONG", "PLAUSIBLE"}


def parse_judge_text(text):
    """Model reply -> list of {pair_id, verdict, tags, key_drug, rationale}."""
    return json.loads(text[text.index("["):text.rindex("]") + 1])


def judge_keep(verdict):
    return verdict in KEEP


def derive_labels(labels, responses):
    """Apply the label rule: curated value wins; otherwise a rare/common (yes/no) model answer
    replaces the rule value for diseases that were sent to the model; 'uncertain' keeps the rule
    value; malignant is 0 when not a neoplasm. Rule values are the shipped values of rows the model
    did not decide."""
    F = labels.merge(responses, on="disease_id", how="left")
    out = []
    for r in F.itertuples(index=False):
        v, rs = r.is_rare, "rule"
        if r.rare_source == "curated":
            rs = "curated"
        elif r.asked_rare and r.response_is_rare in ("rare", "common"):
            v, rs = int(r.response_is_rare == "rare"), "agent"
        c, m, cs = r.is_cancer, r.is_malignant, "rule"
        if r.asked_cancer and r.response_is_cancer in ("yes", "no"):
            c, cs = int(r.response_is_cancer == "yes"), "agent"
            if r.response_is_malignant in ("yes", "no"):
                m = int(r.response_is_malignant == "yes")
            elif not c:
                m = 0
        if not c:
            m = 0
        out.append((r.disease_id, v, c, m, rs, cs))
    return pd.DataFrame(out, columns=["disease_id", "is_rare", "is_cancer", "is_malignant",
                                      "rare_source", "cancer_source"])


def main():
    bad = 0
    J = pd.read_parquet(D / "ixidn_judge_verdicts.parquet")
    n = int((J.verdict.map(judge_keep) != J.keep).sum())
    bad += n
    print(f"judge: {len(J)} pairs, kept {int(J.keep.sum())}, keep mismatches {n}")
    rs = J[J.verdict_source == "expanded_pass_rescored"]
    agree = int((rs.prior_verdict.map(judge_keep) == rs.verdict.map(judge_keep)).sum())
    print(f"re-judged {len(rs)} pairs, keep/drop agreement {agree}/{len(rs)} = {agree / len(rs):.4f}")

    L = pd.read_parquet(D / "disease_labels.parquet")
    A = pd.read_parquet(D / "disease_label_responses.parquet")
    P = derive_labels(L, A).set_index("disease_id").loc[L.disease_id]
    for col in ["is_rare", "is_cancer", "is_malignant", "rare_source", "cancer_source"]:
        n = int((P[col].values != L[col].values).sum())
        bad += n
        print(f"labels: {col} mismatches {n}")
    n = int(((L.is_rare == 1) != (L.disease_class == "rare")).sum())
    bad += n
    print(f"labels: {len(L)} diseases, {int((L.rare_source == 'agent').sum())} rare and "
          f"{int((L.cancer_source == 'agent').sum())} cancer labels set by the model; "
          f"disease_class mismatches {n}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
