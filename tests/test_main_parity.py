"""Parity test: every number in the DisAtlas main text (Abstract, Introduction, Results, Discussion, Methods).

tests/main_text_numbers.tsv lists each number token of the manuscript main text with its kind:
  computed     reproduced here from repo data (data/derived) or from the figure/table value tables written by
               figures/main/*.py and tables/main/*.py; asserted at the printed precision (op eq), as an upper
               bound (op le, "rank N or better") or within a stated tolerance (op approx, Monte-Carlo quantities)
  data_track   network-construction and annotation counts; the channel, graph and benchmark counts among them are
               asserted by tests/test_broad_channel.py, tests/test_specific_channel.py and tests/test_data_model.py
               (and, in COLD mode, tests/test_cold_exact.py)
  si_track     numbers quoted from Supplementary analyses, asserted by the SI parity tests
  parameter, notation, version, date, year, name, identifier, history
               settings, conventions and citations; not data-derived, kept in the registry so every token is
               accounted for
Run: pytest tests/test_main_parity.py   (PARITY_B=2000 bootstrap draws by default; the figures must be built first)
"""
import csv, glob, math, os, re, sys
import pytest
from tolerance import DECIMAL_TOL, dtol_ok

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import main_quantities as MQ

REG = list(csv.DictReader(open(os.path.join(HERE, "main_text_numbers.tsv")), delimiter="\t"))
KINDS = {"computed", "data_track", "si_track", "parameter", "notation", "version", "date", "year", "name", "identifier", "history", "erratum"}


def _value_tables():
    q = {}
    for f in glob.glob(os.path.join(ROOT, "figures", "main", "*_values.tsv")) + glob.glob(os.path.join(ROOT, "tables", "main", "*_values.tsv")):
        nm = os.path.basename(f).replace("_values.tsv", "")
        for r in csv.reader(open(f), delimiter="\t"):
            if len(r) >= 2:
                try: q[f"{nm}:{r[0]}"] = float(r[-1])
                except ValueError: pass
    for f in glob.glob(os.path.join(ROOT, "figures", "main", "fig5_panel_c_*.tsv")):
        nm = os.path.basename(f).replace(".tsv", "")
        for r in csv.DictReader(open(f), delimiter="\t"):
            for c, v in r.items():
                try: q[f"{nm}:{r['stratum']}|{c}"] = float(v)
                except (ValueError, TypeError): pass
    import json
    pats = ["figures/**/*.json", "tables/**/*.json", "data/derived/**/*.json", "data/zenodo/graph/*.json"]
    for f in sorted({g for p in pats for g in glob.glob(os.path.join(ROOT, p), recursive=True)}):
        if os.path.getsize(f) > 5e6:
            continue
        nm = os.path.relpath(f, ROOT)[:-5]

        def walk(k, o):
            if isinstance(o, dict):
                for kk, v in o.items(): walk(f"{k}.{kk}" if k else str(kk), v)
            elif isinstance(o, list) and len(o) <= 4:
                for i, v in enumerate(o): walk(f"{k}[{i}]", v)
            elif isinstance(o, (int, float)) and not isinstance(o, bool):
                q[f"{nm}:{k}"] = float(o)
        try:
            walk("", json.load(open(f)))
        except ValueError:
            pass
    return q


@pytest.fixture(scope="module")
def Q():
    q = MQ.all_quantities(); q.update(_value_tables()); return q


def _parse(tok):
    t = tok.replace(",", "").replace("\u2212", "-").replace("+", "").strip("[]();:. ")
    pct = t.endswith("%"); t = t.rstrip("%")
    dec = len(t.split(".")[1]) if "." in t else 0
    return float(t), dec, pct


def test_registry_complete():
    assert len(REG) > 0 and all(r["kind"] in KINDS for r in REG), {r["kind"] for r in REG} - KINDS
    assert len({r["id"] for r in REG}) == len(REG)
    assert all(r["key"] for r in REG if r["kind"] == "computed")
def test_no_stale_pending():
    """a PENDING row must become a plain computed row once its value table provides the key"""
    q = _value_tables()
    stale = [r["id"] for r in REG if r["note"].startswith("PENDING") and r["key"] in q]
    assert not stale, stale


COMPUTED = [r for r in REG if r["kind"] == "computed"]


@pytest.mark.parametrize("row", COMPUTED, ids=[f"{r['id']}:{r['token']}" for r in COMPUTED])
def test_main_text_number(row, Q):
    if row["note"].startswith("PENDING") and row["key"] not in Q:
        pytest.skip(row["note"])  # value table not yet regenerated for the current manuscript
    assert row["key"] in Q, f"quantity {row['key']} not computed"
    v = Q[row["key"]]; t, dec, pct = _parse(row["token"])
    if pct and abs(v) <= 1 and abs(t) > 1: v = 100 * v
    if row["op"] == "le":
        assert t >= round(v, dec), (row["context"], v)
    elif row["op"] == "ge":  # 'at or above <token>': the printed bound may not exceed the value
        assert t <= v + 1e-12, (row["context"], v)
    elif row["op"] == "dtol":
        assert dtol_ok(v, t, row["token"]), (row["context"], v, DECIMAL_TOL)
    elif row["op"] == "approx":
        assert abs(v - t) <= float(row["tol"]), (row["context"], v)
    else:
        assert round(v, dec) == pytest.approx(t, abs=1e-9) or abs(v - t) <= 0.5 * 10 ** -dec + 1e-9, (row["context"], v)


ERRATA = [r for r in REG if r["kind"] == "erratum"]


def _match(v, tok):
    t, dec, pct = _parse(tok)
    if pct and abs(v) <= 1 and abs(t) > 1: v = 100 * v
    return round(v, dec) == pytest.approx(t, abs=1e-9) or abs(v - t) <= 0.5 * 10 ** -dec + 1e-9


@pytest.mark.parametrize("row", ERRATA, ids=[f"{r['id']}:{r['token']}" for r in ERRATA])
def test_main_text_erratum(row, Q):
    """A number of the final manuscript that docs/ERRATA.md lists as stale. `key` is the corrected quantity from the
    released outputs; the printed token must differ from it (once the text is corrected the row becomes `computed`).
    When the note names `stale=<key>`, the printed token is also reproduced from that shipped stale source."""
    assert row["note"].startswith("ERRATUM E"), row["note"]
    assert row["key"] in Q, f"quantity {row['key']} not computed"
    assert not _match(Q[row["key"]], row["token"]), "printed value now matches the corrected quantity: make the row computed"
    m = re.search(r"stale=(\S+)", row["note"])
    if m:
        assert m.group(1) in Q and _match(Q[m.group(1)], row["token"]), (m.group(1), Q.get(m.group(1)))
