"""Shared helpers for the table scripts under tables/main/ and tables/SI/.

Every table script computes each cell from repo data, writes <name>.csv and <name>.md
next to itself, and checks every cell against the value printed in the manuscript
(EXPECTED in the script). A mismatch raises, so a passing run is a parity check.
"""
from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
RAW = REPO_ROOT / "data" / "00_raw"
DATA = REPO_ROOT / "data"
INPUTS = REPO_ROOT / "tables" / "inputs"


def f3(x):
    return f"{x:.3f}"


def signed3(x):
    s = f"{x:+.3f}"
    return "+0.000" if s == "-0.000" else s


def thousands(n):
    return f"{int(n):,}"


def pct1(x):
    return f"{100 * x:.1f}%"


def to_markdown(df):
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in df.iterrows():
        out.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    return "\n".join(out) + "\n"


def write_and_check(df, script_file, name, expected=None, caption="", known=()):
    """Write CSV + Markdown beside the calling script and assert parity with EXPECTED.

    `known` lists (row, column) cells where the printed manuscript value is a documented
    manuscript error (see tables/README.md); they are reported, not silently accepted."""
    here = Path(script_file).resolve().parent
    df.to_csv(here / f"{name}.csv", index=False)
    md = to_markdown(df)
    if caption:
        md += "\n" + caption.strip() + "\n"
    (here / f"{name}.md").write_text(md)
    if expected is not None:
        exp = pd.DataFrame(expected, columns=df.columns).astype(str)
        got = df.astype(str).reset_index(drop=True)
        bad = [(i, c, exp.at[i, c], got.at[i, c]) for i in range(len(exp)) for c in exp.columns
               if exp.at[i, c] != got.at[i, c]]
        for b in [b for b in bad if (b[0], b[1]) in set(known)]:
            print("KNOWN MANUSCRIPT DISCREPANCY row %d col %r: manuscript=%r rebuilt=%r" % b)
        bad = [b for b in bad if (b[0], b[1]) not in set(known)]
        if len(exp) != len(got) or bad:
            for b in bad:
                print("MISMATCH row %d col %r: manuscript=%r rebuilt=%r" % b)
            raise SystemExit(f"{name}: parity FAILED")
        tail = f" except {len(known)} documented manuscript discrepancies" if known else ""
        print(f"{name}: {len(got)} rows x {len(got.columns)} cols, all cells equal the manuscript{tail}")
    return df
