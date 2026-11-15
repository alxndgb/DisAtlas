"""Every tracked data file is read by shipped code, and data/DATA_INDEX.tsv plus the data/README.md "Read by" block are current
(regenerate with: python scripts/maintainer/data_index.py)."""
import importlib.util
from pathlib import Path

spec = importlib.util.spec_from_file_location("data_index", Path(__file__).resolve().parents[1] / "scripts/maintainer/data_index.py")
DI = importlib.util.module_from_spec(spec); spec.loader.exec_module(DI)


def test_index_current_and_no_dangling_data():
    stale, unread = DI.main(check=True)
    assert not unread, f"data files read by no shipped code: {unread}"
    assert not stale, f"stale: {stale}; run python scripts/maintainer/data_index.py"
