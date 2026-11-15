"""scripts/fetch_open_targets.py (optional Figure 1 fetcher): import and --help only, no network access."""
import importlib.util
import socket
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "fetch_open_targets.py"


@pytest.fixture
def no_network(monkeypatch):
    def refuse(*a, **k):
        raise AssertionError("network access attempted")
    monkeypatch.setattr(socket, "socket", refuse)
    monkeypatch.setattr(socket, "create_connection", refuse)


def load():
    spec = importlib.util.spec_from_file_location("fetch_open_targets", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_import_has_public_api_only(no_network):
    mod = load()
    assert mod.API_URL == "https://api.platform.opentargets.org/api/v4/graphql"
    assert callable(mod.main)


@pytest.mark.parametrize("argv", [["--help"], ["fetch", "--help"], ["list-diseases", "--help"], ["merge", "--help"]])
def test_help(no_network, monkeypatch, tmp_path, capsys, argv):
    monkeypatch.setenv("OT_DATA_ROOT", str(tmp_path))
    mod = load()
    with pytest.raises(SystemExit) as e:
        mod.main(argv)
    assert e.value.code == 0
    assert "usage" in capsys.readouterr().out


def test_not_in_pipeline():
    assert "fetch_open_targets" not in (SCRIPT.parents[1] / "run_pipeline.sh").read_text()
