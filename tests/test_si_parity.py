"""Every number of the final Supplementary Information sections that no other parity test asserts is registered in
tests/si_text_numbers.tsv (same schema as tests/main_text_numbers.tsv plus `classified_by`).

kind:
  computed   `key` is a quantity of tests/si_quantities.py; the printed token must equal it at its printed precision
             (op eq), or satisfy op le / ge / approx (tol).
  erratum    docs/ERRATA.md lists the printed token as stale; `key` is the regenerated value, which must differ from the
             token; `stale=<key>` in the note, when present, must reproduce the printed token.
  pending    the source is not shipped yet; the note names what is missing (skipped, never passed silently).
  parameter, notation, label, citation, year, design: not checked against data; `note` gives the reason.
classified_by records which pass proposed the label and whether the session model reviewed or changed it.
"""
import csv, os, re
import pytest
from tolerance import DECIMAL_TOL, dtol_ok
import si_quantities as SQ

HERE = os.path.dirname(os.path.abspath(__file__))
REG = list(csv.DictReader(open(os.path.join(HERE, "si_text_numbers.tsv")), delimiter="\t", quoting=csv.QUOTE_NONE))
KINDS = {"computed", "erratum", "pending", "parameter", "notation", "label", "citation", "year", "design"}
COLS = ["id", "section", "line", "token", "kind", "key", "op", "tol", "note", "context", "classified_by"]


@pytest.fixture(scope="module")
def Q():
    return SQ.all_si_quantities_with_derived()


def _parse(tok):
    t = tok.replace(",", "").replace("\u2212", "-").replace("+", "").strip("[]();:. ")
    pct = t.endswith("%"); t = t.rstrip("%")
    dec = len(t.split(".")[1]) if "." in t else 0
    return float(t), dec, pct


def _match(v, tok):
    t, dec, pct = _parse(tok)
    if pct and abs(v) <= 1 and abs(t) > 1: v = 100 * v
    return round(v, dec) == pytest.approx(t, abs=1e-9) or abs(v - t) <= 0.5 * 10 ** -dec + 1e-9


def test_registry_complete():
    assert REG and list(REG[0].keys()) == COLS
    assert all(r["kind"] in KINDS for r in REG), {r["kind"] for r in REG} - KINDS
    assert len({r["id"] for r in REG}) == len(REG)
    assert all(r["key"] for r in REG if r["kind"] in ("computed", "erratum"))
    assert all(r["note"] for r in REG if r["kind"] in ("pending", "erratum", "design", "citation"))
    assert all(r["classified_by"] for r in REG)


def test_no_stale_pending(Q):
    stale = [r["id"] for r in REG if r["kind"] == "pending" and r["key"] and r["key"] in Q]
    assert not stale, stale


COMPUTED = [r for r in REG if r["kind"] == "computed"]


@pytest.mark.parametrize("row", COMPUTED, ids=[f"{r['id']}:{r['token']}" for r in COMPUTED])
def test_si_text_number(row, Q):
    assert row["key"] in Q, f"quantity {row['key']} not computed"
    v = Q[row["key"]]; t, dec, pct = _parse(row["token"])
    if pct and abs(v) <= 1 and abs(t) > 1: v = 100 * v
    if row["op"] == "le":
        assert t >= round(v, dec), (row["context"], v)
    elif row["op"] == "ge":
        assert t <= v + 1e-12, (row["context"], v)
    elif row["op"] == "dtol":
        assert dtol_ok(v, t, row["token"]), (row["context"], v, DECIMAL_TOL)
    elif row["op"] == "approx":
        assert abs(v - t) <= float(row["tol"]), (row["context"], v)
    else:
        assert _match(v, row["token"]), (row["context"], v)


ERRATA = [r for r in REG if r["kind"] == "erratum"]


@pytest.mark.parametrize("row", ERRATA, ids=[f"{r['id']}:{r['token']}" for r in ERRATA])
def test_si_text_erratum(row, Q):
    assert row["note"].startswith("ERRATUM E"), row["note"]
    assert row["key"] in Q, f"quantity {row['key']} not computed"
    assert not _match(Q[row["key"]], row["token"]), "printed value now matches: make the row computed"
    m = re.search(r"stale=(\S+)", row["note"])
    if m:
        assert m.group(1) in Q and _match(Q[m.group(1)], row["token"]), (m.group(1), Q.get(m.group(1)))


PENDING = [r for r in REG if r["kind"] == "pending"]


@pytest.mark.parametrize("row", PENDING, ids=[f"{r['id']}:{r['token']}" for r in PENDING])
def test_si_text_pending(row):
    pytest.skip(row["note"])
