"""SI event-time odds-ratio tables and prose numbers: regenerate and compare with the committed tables and manuscript sentences."""
import json, subprocess, sys
from pathlib import Path
import pytest
ROOT = Path(__file__).resolve().parents[1]
SPEC = json.loads((ROOT / "tests" / "si_event_time_or.json").read_text(encoding="utf-8"))

@pytest.fixture(scope="module")
def built(tmp_path_factory):
    out = tmp_path_factory.mktemp("si_or")
    subprocess.run([sys.executable, str(ROOT / "figures" / "SI" / "si_event_time_or_tables.py"), str(out)], check=True, cwd=ROOT)
    return out

def test_tables_byte_identical(built):
    assert (built / "SI_event_time_or_tables.md").read_bytes() == (ROOT / "figures" / "SI" / "SI_event_time_or_tables.md").read_bytes()

@pytest.mark.parametrize("key", sorted(SPEC["sentences"]))
def test_prose_numbers(built, key):
    P = json.loads((built / "SI_event_time_or_tables_numbers.json").read_text())["prose"]
    n21 = {**SPEC.get("final_model_note02", {}), **SPEC.get("text_defects", {})}
    if key in n21:
        # printed value differs from the final model or the printed sentence is defective: pin the final value
        assert P[key] == n21[key]["final"], (key, P[key], n21[key]["final"])
        return
    sent = SPEC["sentences"][key]
    for s in P[key]: assert s in sent, (key, s, sent)


def test_note02_entries_documented():
    for k, v in {**SPEC.get("final_model_note02", {}), **SPEC.get("text_defects", {})}.items():
        assert v["printed"] and v["final"] and v["reason"], k
