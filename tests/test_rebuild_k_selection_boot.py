"""Figure S16 a: bootstrap replicate 0 of the K-selection table is rebuilt at every K from the released cosine matrix
(disatlas.rebuild.k_selection_boot). Needs the optional cluster extra (requirements-cluster.txt); skipped otherwise.
The full 100-replicate rebuild is  python -m disatlas.rebuild.k_selection_boot OUTDIR  (hours)."""
from pathlib import Path
import numpy as np, pandas as pd
import pytest

pytest.importorskip("kmedoids")
from disatlas.rebuild import k_selection_boot as KS

STORED = Path(__file__).resolve().parents[1] / "data/derived/si_cluster/k_selection_boot.tsv.gz"


def test_replicate0():
    S = pd.read_csv(STORED, sep="\t"); S = S[S.rep == 0].reset_index(drop=True)
    T = KS.compute(reps=[0])
    assert list(T.columns) == list(S.columns)
    for c in ("kind", "rep", "K", "min_size", "sub_seed"):
        assert (T[c].to_numpy() == S[c].to_numpy()).all(), c
    for c in ("silhouette", "ari_vs_full"):
        assert np.abs(T[c].to_numpy(float) - S[c].to_numpy(float)).max() <= 1e-9, (c, np.abs(T[c] - S[c]).max())


@pytest.mark.skipif(not __import__("os").environ.get("DISATLAS_KSEL_RUN"), reason="set DISATLAS_KSEL_RUN to a sharded run dir")
def test_full_100_replicates():
    """All 100 replicates x 15 K from scripts/slurm/k_selection_boot_array.sh (150 tasks), bit-equal to the shipped table."""
    S = pd.read_csv(STORED, sep="\t").sort_values(["rep", "K"]).reset_index(drop=True)
    T = KS.merge(__import__("os").environ["DISATLAS_KSEL_RUN"])
    assert T.shape == S.shape == (1500, 7)
    for c in S.columns:
        assert (T[c].to_numpy() == S[c].to_numpy()).all(), c
