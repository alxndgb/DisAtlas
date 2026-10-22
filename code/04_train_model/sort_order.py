"""Original disease-ID sort order for blinded IDs.

Shipped disease IDs for non-public vocabularies are blinded (DIS_B...), so sorting them as strings no
longer gives the order the frozen model was trained with. data/03_rebuilt/frozen/disease_sort_rank.parquet
stores each disease's rank in the original string order; sorting by it reproduces that order exactly.
"""
import os

import pandas as pd

_RANK_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "data", "03_rebuilt", "frozen",
                          "disease_sort_rank.parquet")
_RANK = None


def _rank():
    global _RANK
    if _RANK is None:
        r = pd.read_parquet(_RANK_PATH)
        _RANK = dict(zip(r["disease_id"], r["sort_rank"]))
    return _RANK


def original_order(ids):
    """Sort disease IDs in the order used to train the frozen model."""
    rk = _rank()
    missing = [i for i in ids if i not in rk]
    if missing:
        raise KeyError(f"{len(missing)} disease IDs have no stored sort rank, e.g. {missing[:3]}")
    return sorted(ids, key=rk.__getitem__)
