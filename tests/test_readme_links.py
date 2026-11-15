"""Every relative link in README.md resolves to a tracked file or directory (fast; no data needed)."""
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
README = ROOT / "README.md"

LINK = re.compile(r"(?<!!)\[[^\]\n]*\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)")
HREF = re.compile(r'<a\s+href="([^"]+)"')
FENCE = re.compile(r"^(```|~~~)")


def _readme_links():
    links, fenced = [], False
    for n, line in enumerate(README.read_text(encoding="utf-8").splitlines(), 1):
        if FENCE.match(line.strip()):
            fenced = not fenced
            continue
        if fenced:
            continue
        line = re.sub(r"`[^`]*`", "", line)  # ignore inline code
        for m in LINK.finditer(line):
            links.append((n, m.group(1)))
        for m in HREF.finditer(line):  # the item maps are HTML tables (column widths), links are <a href>
            links.append((n, m.group(1)))
    return links


def _tracked():
    out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:  # not a git checkout (e.g. a release tarball): fall back to existence
        return None
    return set(filter(None, out.stdout.split("\0")))


def test_readme_relative_links_resolve_to_tracked_files():
    tracked = _tracked()
    links = _readme_links()
    assert links, "README has no links?"
    bad = []
    for n, target in links:
        if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.I) or target.startswith(("#", "<")):
            continue  # external URL, mailto, in-page anchor
        path = target.split("#", 1)[0].split("?", 1)[0]
        if not path:
            continue
        rel = (README.parent / path).resolve().relative_to(ROOT).as_posix()
        if not (ROOT / rel).exists():
            bad.append(f"line {n}: {target} (missing)")
        elif tracked is not None:
            is_dir = (ROOT / rel).is_dir()
            ok = any(t.startswith(rel.rstrip("/") + "/") for t in tracked) if is_dir else rel in tracked
            if not ok:
                bad.append(f"line {n}: {target} (not tracked in git)")
    assert not bad, "broken README links:\n" + "\n".join(bad)


def test_readme_item_tables():
    # The three item maps are HTML tables (GitHub keeps width attributes on th, not style): columns are
    # Item | Content | Script | Output; the item label and its png/pdf/md/tsv links share the first column,
    # Content is the widest, every Figure/Table row has the label followed by at least one bracketed link
    # in column 1 ('Figure 1 [png] [pdf]', no separator),
    # and Script/Output cells show file names only (the full path is the link target).
    text = README.read_text(encoding="utf-8")
    tables = re.findall(r"<table>.*?</table>", text, re.S)
    assert len(tables) == 3, len(tables)  # main text, supplementary figures, supplementary tables
    nrows, missing, long_paths, bad_first = 0, [], [], []
    for tb in tables:
        ths = re.findall(r'<th width="(\d+)%">([^<]+)</th>', tb)
        names = [n for _, n in ths]
        assert names == ["Item", "Content", "Script", "Output"], names
        w = {n: int(x) for x, n in ths}
        assert sum(w.values()) == 100, w
        assert w["Content"] >= 40 and w["Content"] > max(w["Item"], w["Script"], w["Output"]), w
        for tr in re.findall(r"<tr>(.*?)</tr>", tb):
            cells = re.findall(r"<td>(.*?)</td>", tr)
            if not cells:
                continue
            assert len(cells) == 4, cells
            first = cells[0]
            item = first.split(" <a ", 1)[0]
            if "·" in first or "<" in item:
                bad_first.append(first)  # label, a space, then '<a href="...">[png]</a> <a ...>[pdf]</a>', no separator
            links = first[len(item):]
            if links and not re.fullmatch(r'( <a href="[^"]+"( title="[^"]*")?>\[(png|pdf|md|tsv|csv|parquet|json)\]</a>)+', links):
                bad_first.append(first)
            if re.match(r"(Figure|Table|Tables) ", item):
                nrows += 1
                if not links:
                    missing.append(item)
            for c in cells[2:]:
                for word in re.sub(r"<[^>]+>", "", c).split():
                    if re.match(r"\(?(figures|tables|data|tests|src|scripts|docs|release|notebooks|\.\.)/", word) and "*" not in word:
                        long_paths.append((item, word))
    assert nrows >= 40, nrows
    assert not missing, f"item rows without label + link in column 1: {missing}"
    assert not bad_first, f"first cell must be 'Figure 1 [png] [pdf]' style: {bad_first[:3]}"
    assert not long_paths, f"path text in Script/Output cells (show file names only): {long_paths}"
