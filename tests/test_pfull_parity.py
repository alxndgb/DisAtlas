"""P_full parity record: replaying the authors' P_full retrain (shipped per-seed cosines of the IxIDN-Matched
pairs, data/zenodo/model_outputs/pfull_reference) against the released model gives the recorded verdict."""
import json
import os
import sys

import pytest

from disatlas import paths

sys.path.insert(0, str(paths.ROOT / "scripts"))
import pfull_parity as PP  # noqa: E402

needs = pytest.mark.skipif(not (os.path.exists(PP.REC) and os.path.isdir(PP.RELEASED)), reason="Zenodo bundle not present")


@needs
def test_pfull_parity_replay():
    r = PP.parity(PP.table())
    rec = json.load(open(PP.RECORD))
    assert r["verdict"] == rec["verdict"] == "PARITY"
    for k, v in rec.items():
        if isinstance(v, float):
            assert r[k] == pytest.approx(v, abs=1e-12), k
        else:
            assert r[k] == v, k
    assert r["n"] == r["shared"] == 1094 and r["ci_lo"] > PP.MARGIN
