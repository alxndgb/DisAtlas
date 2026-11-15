"""Overlay of corrected values on the frozen manuscript target.

`release/manuscript_overlay.tsv` lists every printed value that the correction notes 01-04 replace
(label, location, kind, old value, new value, note). Table rows are the reproduction target: each new
row must appear in the regenerated Supplementary table, and no old row may remain there.
"""
import csv, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OVL = ROOT / "release" / "manuscript_overlay.tsv"
TABLES = {
    "tbl:SI_event_time_deltas": "figures/SI/SI_event_time_tables.md",
    "tbl:SI_event_time_subsets": "figures/SI/SI_event_time_tables.md",
    "tbl:SI_event_time_or_yearly": "figures/SI/SI_event_time_or_tables.md",
    "tbl:SI_event_time_auc": "figures/SI/SI_event_time_or_tables.md",
    "tbl:SI_event_time_auc_yearly": "figures/SI/SI_event_time_or_tables.md",
    "tbl:SI_fda_dated_pairs": "figures/SI/SI_fda_prioritisation_text.md",
}


def nrm(s):
    s = re.sub(r'\{custom-style="[^"]*"\}', "", s)
    return re.sub(r"[\[\]\*\s|:\-\u2013\u2212]", "", s)


def rows():
    with open(OVL) as fh:
        return list(csv.DictReader((l for l in fh if not l.startswith("#")), delimiter="\t"))


def test_overlay_header_and_columns():
    head = OVL.read_text().split("\n")[:2]
    assert "reproduction target" in head[1]
    r = rows()
    assert r and set(r[0]) == {"label", "location", "kind", "old_value", "new_value", "note"}
    assert all(x["old_value"] != x["new_value"] for x in r)


def test_overlay_table_rows_regenerated():
    r = [x for x in rows() if x["kind"] == "table row"]
    assert {x["label"] for x in r} == set(TABLES)
    for lab, f in TABLES.items():
        got = {nrm(l) for l in (ROOT / f).read_text().split("\n") if l.strip().startswith("|")}
        sel = [x for x in r if x["label"] == lab]
        miss = [x["new_value"][:80] for x in sel if nrm(x["new_value"]) not in got]
        stale = [x["old_value"][:80] for x in sel if nrm(x["old_value"]) in got]
        assert not miss, (lab, len(miss), miss[:3])
        assert not stale, (lab, len(stale), stale[:3])
