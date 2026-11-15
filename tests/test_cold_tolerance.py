"""Cold mode: the ensemble retrained on the P_full graph (build/model) passes the pre-registered parity test
against the released model (scripts/pfull_parity.py: held-out IxIDN-Matched test pairs, cosine readout,
5 seeds; disease-cluster bootstrap 95 % CI lower bound of the AUROC difference > -0.02 and coverage >= 95 %).
The released model was trained on a graph whose specific channel came from BIKG inputs that are not released,
so a retrain is judged by this parity rule, not by score identity."""
import sys

import pytest

from disatlas import paths

pytestmark = pytest.mark.skipif(paths.MODE != "cold", reason="cold mode only (DISATLAS_DATA_MODE=cold)")
M = paths.BUILD / "model" / paths.MODEL_NAME


def test_cold_retrain_parity():
    if not all((M / "embeddings" / f"node2vec_s{s}.npz").exists() for s in (4707, 4708, 4709, 4710, 4711)):
        pytest.skip(f"{M} not built (retrain stage not run)")
    sys.path.insert(0, str(paths.ROOT / "scripts"))
    import pfull_parity as PP
    r = PP.parity(PP.table(model=str(M)))
    assert r["coverage"] >= PP.COVER and r["ci_lo"] > PP.MARGIN, r
