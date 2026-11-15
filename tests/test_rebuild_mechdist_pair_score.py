"""Parity: disatlas.rebuild.mechdist_pair_score vs data/derived/si_mechdist_pair_score (Supplementary Figure S6).
Class (c) for one value: the shared-neighbour count (52), the no-direct-edge flag, the hub membership, the Youden threshold,
test sensitivity / specificity, the cosine table and the six schematic cosines reproduce; the hub degree of the 'neoplasm'
umbrella term is 3,006 in the staged training graph against 3,036 shipped (shipped value kept)."""
import json
import numpy as np, pandas as pd
from disatlas.rebuild import mechdist_pair_score as M

O = json.load(open(M.REF / "s6_inputs.json"))


def test_cosine_table():
    C = M.cosines(); R = pd.read_csv(M.REF / "ixidn_matched_cosine_all.tsv", sep="\t")
    assert (C.y.to_numpy() == R.y.to_numpy()).all() and (C.pool.to_numpy() == R.pool.to_numpy()).all()
    np.testing.assert_allclose(C.c.to_numpy(), R.c.to_numpy(), rtol=0, atol=2e-6)


def test_schematic_cosines():
    R = pd.read_csv(M.REF / "schematic_cd_control_cosines.tsv", sep="\t")
    assert np.allclose(M.schematic(R).cosine_score.to_numpy(), R.cosine_score.to_numpy(), rtol=0, atol=1e-9)


def test_threshold_and_test_operating_point():
    thr, se, sp = M.threshold_sens_spec(M.cosines())
    assert round(thr, 6) == O["threshold_youden_tune"]
    assert abs(se - O["test_sensitivity"]) < 1e-12 and abs(sp - O["test_specificity"]) < 1e-12


def test_training_graph_counts():
    T, sh = M.graph_counts(O["pair"]); R = O["training_graph"]
    assert T["n_shared"] == R["n_shared"] and T["direct_AB"] == R["direct_AB"] and T["hub_in_shared"] == R["hub_in_shared"]
    assert set(O["figure_shared"]) <= sh
    assert T["hub_deg"] == 3006 and R["hub_deg"] == 3036      # (c): shipped hub degree not reproduced, see module docstring
