"""Comparator scores: each `system|run` column of data/raw/benchmarks/comparators/ixidn_matched_comparator_scores.tsv
is the DDA column of the stored per-run comparator file (exact match on blind_id = pair_id, checked once against the
originals on 2026-11-14; the originals are not shipped). This test pins each column by checksum and checks the set of
columns against the systems and runs the manuscript cites."""
import hashlib, json
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
F = ROOT / "data/raw/benchmarks/comparators/ixidn_matched_comparator_scores.tsv"
REF = json.loads((Path(__file__).parent / "parity/comparator_column_sha256.json").read_text())["columns"]
CITED = {"Claude-Sci Opus 5.5": 5, "Claude-Chat Opus 5.5": 5, "Gemini Pro 3.1": 3, "TxGNN published ckpt": 1, "Phylo": 5}


def _m():
    return pd.read_csv(F, sep="\t", dtype={"pair_id": str})


def _sha(v):
    return hashlib.sha256("\n".join("NA" if pd.isna(x) else repr(float(x)) for x in v).encode()).hexdigest()


def test_comparator_columns_checksum():
    M = _m().sort_values("pair_id")
    assert len(M) == 4364 and M.pair_id.is_unique
    assert [c for c in M.columns if c != "pair_id"] == list(REF)
    bad = [c for c, h in REF.items() if _sha(M[c].values) != h]
    assert not bad, bad


def test_comparator_runs_cited():
    cols = [c for c in _m().columns if c != "pair_id"]
    n = pd.Series([c.split("|")[0] for c in cols]).value_counts().to_dict()
    assert n == CITED, n
    assert all("Deepwalk" in c for c in cols if c.startswith("Phylo|"))   # no unused Phylo variants
    assert not any("argon" in c.lower() for c in cols)                    # no Gemini-Argon runs


def test_txgnn_unmapped_ranked_last():
    M = _m().set_index("pair_id"); tx = M["TxGNN published ckpt|DistMult"]
    assert int(tx.isna().sum()) == 213                                     # unmappable pairs are NA in the raw scores
    F5 = pd.read_csv(ROOT / "data/derived/fig5_clustering/fig5c_txgnn_pair_scores.tsv", sep="\t",
                     dtype={"blind_id": str}).set_index("blind_id").score.reindex(tx.index)
    m = tx.notna()
    assert (F5[m] == tx[m]).all()
    assert (F5[~m] < F5[m].min()).all() and F5[~m].nunique() == 1           # one sentinel below every scored pair
