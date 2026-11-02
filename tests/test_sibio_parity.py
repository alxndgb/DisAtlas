"""Parity tests for the SI biology items.

Each SI biology item has a spec file tests/parity/sibio_<group>.json:
    {"items": [{"label": "SI_...", "script": "figures/SI/si_x.py",
                "numbers": "figures/SI/si_x_numbers.json",
                "status": "exact" | "visual-only" | "precomputed-constant",
                "reason": "...",                          # required unless status == exact
                "checks": [{"key": "a.b.c", "expected": 0.123, "tol": 5e-4, "ms": "quoted manuscript text"},
                           {"key": "a.p", "op": "lt", "bound": 1e-30, "ms": "p < 10^-30"}],
                "text_checks": [{"file": "figures/SI/x.md", "numbers": ["1566", "5.73"], "ms": "..."}]}],
     "unchecked": [{"token": "...", "reason": "...", "ms": "..."}],          # not numerical results (labels, ordinals)
     "manuscript_mismatch": [{"numbers": "...json", "key": "...", "data_value": ..., "ms_value": "...",
                              "reason": "..."}]}                            # printed value differs from the regenerated one
``key`` is a dotted path into the numbers json (list indices allowed). ``expected`` is the value as printed in
the manuscript; ``tol`` is half a unit of the last printed digit unless stated otherwise. Strings and integers
compare exactly. ``text_checks`` require the printed numbers to appear, in order and contiguously, among the numbers
of the generated markdown table. ``manuscript_mismatch`` entries pin the regenerated value (so the pipeline cannot
drift) and record the differing printed value for correction in the manuscript. Run the scripts first.
"""
import json
import operator
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SPECS = sorted((REPO / "tests" / "parity").glob("sibio_*.json"))
ITEMS = [(s.stem, it) for s in SPECS for it in json.loads(s.read_text())["items"]]
MISMATCH = [(s.stem, m) for s in SPECS for m in json.loads(s.read_text()).get("manuscript_mismatch", [])]
OPS = {"lt": operator.lt, "le": operator.le, "gt": operator.gt, "ge": operator.ge}
NUM = re.compile(r"[−-]?\d[\d,]*(?:\.\d+)?")


def _get(d, key):
    for k in key.split("."):
        d = d[int(k)] if isinstance(d, list) else d[k]
    return d


def _norm(t):
    return t.replace("−", "-").replace(",", "").rstrip("%").lstrip("+")


def _md_numbers(path):
    out = []
    for m in NUM.finditer(path.read_text(encoding="utf-8")):
        a = m.start()
        if a > 0 and (m.string[a - 1].isalpha() or m.string[a - 1] in "_/"):
            continue
        out.append(_norm(m.group(0)))
    return out


def test_specs_exist():
    assert SPECS, "no tests/parity/sibio_*.json spec files"


@pytest.mark.parametrize("group,item", ITEMS, ids=[f"{g}:{i['label']}" for g, i in ITEMS])
def test_item(group, item):
    assert item["status"] in {"exact", "visual-only", "precomputed-constant"}
    if item["status"] != "exact":
        assert item.get("reason"), f"{item['label']}: non-exact status needs a reason"
    assert (REPO / item["script"]).exists(), item["script"]
    if not item.get("checks") and not item.get("text_checks"):
        assert item["status"] == "visual-only", f"{item['label']}: no checks but status {item['status']}"
        return
    bad = []
    if item.get("checks"):
        nums = json.loads((REPO / item["numbers"]).read_text())
        for c in item["checks"]:
            got = _get(nums, c["key"])
            if "op" in c:
                ok = OPS[c["op"]](float(got), float(c["bound"]))
                exp = f"{c['op']} {c['bound']}"
            else:
                exp = c["expected"]
                if isinstance(exp, str) or (isinstance(exp, int) and "tol" not in c):
                    ok = got == exp
                else:
                    ok = abs(float(got) - float(exp)) <= float(c.get("tol", 0)) + 1e-12
            if not ok:
                bad.append(f"{c['key']}: got {got!r}, manuscript {exp!r} ({c.get('ms', '')})")
    cache = {}
    for t in item.get("text_checks", []):
        f = REPO / t["file"]
        nums = cache.setdefault(f, _md_numbers(f))
        seq = [_norm(x) for x in t["numbers"]]
        if not any(nums[i:i + len(seq)] == seq for i in range(len(nums) - len(seq) + 1)):
            bad.append(f"{t['file']}: sequence {t['numbers']} not found ({t.get('ms', '')[:120]})")
    assert not bad, f"{item['label']}:\n" + "\n".join(bad)


@pytest.mark.parametrize("group,m", MISMATCH, ids=[f"{g}:{m['key']}" for g, m in MISMATCH])
def test_manuscript_mismatch_pinned(group, m):
    assert m.get("reason") and "ms_value" in m
    got = _get(json.loads((REPO / m["numbers"]).read_text()), m["key"])
    assert got == m["data_value"], f"{m['key']}: regenerated {got!r} != pinned {m['data_value']!r}"
