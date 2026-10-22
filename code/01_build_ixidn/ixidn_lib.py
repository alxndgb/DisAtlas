"""Shared paths, knobs and helpers for the IxIDN build (raw tables -> DD network -> IxIDN+ / IxIDN-)."""
import json
import os
from collections import defaultdict

import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
RAW = os.path.join(REPO_ROOT, "data", "00_raw")
ZRAW = os.path.join(REPO_ROOT, "data", "zenodo", "00_raw")
OUT = os.path.join(REPO_ROOT, "data", "01_rebuilt")
MONDO_OBO = os.path.join(ZRAW, "ontologies", "mondo.obo")

RARE_SOURCE = "SRC_0B0DBCC4E9"          # rare-disease registry source (layer-5 bypass)
RARE_SUBSETS = ("orphanet_rare", "gard_rare")

FROZEN = os.path.join(RAW, "frozen")
# frozen matching inputs: single constants so a file rename is one edit
MATCH_FEATURES = os.path.join(FROZEN, "ixidn_matched_universe_features.parquet")
MATCH_CANDIDATES = os.path.join(FROZEN, "ixidn_matched_positive_candidates.parquet")
EVAL_UNIVERSE = os.path.join(FROZEN, "ixidn_eval_pair_universe_frozen.tsv.gz")
DISEASE_SPLIT = os.path.join(FROZEN, "disease_split.tsv.gz")
CAUSAL_GENES = os.path.join(FROZEN, "disease_causal_genes.parquet")
JUDGE_VERDICTS = os.path.join(RAW, "llm_annotations", "ixidn_judge_verdicts.parquet")

KNOBS = {
    "drop_decisions": ["exclude"],      # layer 1: non-MoA drugs
    "max_drugs_per_trial": 2,           # layer 0: trial-size cap, counting kept drugs only
    "drug_degree_cap_pct": 99.9,        # layer 3
    "top_k": 40, "max_partners": 100,   # layer 4
    "min_shared_drugs": 2,
    "rare_rare_min_shared": 1, "rare_top_k": 100, "rare_max_partners": 300,
    "coinv_hier_max_children": 10,      # co-investigation veto via direct MONDO parent/child
    "neg_d_min": 3,                     # IxIDN-: pairs >= 3 hops apart in IxIDN+
    "neg_degree_bins": 5,
}


def out_path(name):
    os.makedirs(OUT, exist_ok=True)
    return os.path.join(OUT, name)


def load_triples():
    return pd.read_parquet(os.path.join(RAW, "trials", "trial_triples.parquet"))


def load_diseases():
    return pd.read_parquet(os.path.join(RAW, "diseases", "diseases.parquet"))


def load_drug_decisions():
    return pd.read_parquet(os.path.join(RAW, "drugs", "drug_moa_decision.parquet"))


def mondo_is_a(path=MONDO_OBO):
    """{MONDO term -> set(direct is_a MONDO parents)}."""
    par, cur = defaultdict(set), None
    with open(path) as fh:
        for line in fh:
            if line.startswith("[Term]"):
                cur = None
            elif line.startswith("[Typedef]"):
                cur = "TYPEDEF"
            elif line.startswith("id: ") and cur != "TYPEDEF":
                cur = line[4:].strip()
            elif line.startswith("is_a: ") and cur and cur.startswith("MONDO:"):
                p = line.split()[1]
                if p.startswith("MONDO:"):
                    par[cur].add(p)
    return par


def mondo_rare(path=MONDO_OBO):
    """Set of MONDO ids carrying a rare subset (orphanet_rare | gard_rare)."""
    out, cid = set(), None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("[Term]"):
                cid = None
            elif line.startswith("id: MONDO:"):
                cid = line[4:].strip()
            elif cid and line.startswith("subset: ") and line.split()[1] in RARE_SUBSETS:
                out.add(cid)
    return out


def read_out(name):
    return pd.read_parquet(out_path(name))
