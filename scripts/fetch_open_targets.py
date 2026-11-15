#!/usr/bin/env python3
"""Cached, resumable fetcher for Open Targets Platform disease-target associations.

Pages ``disease(efoId).associatedTargets`` from the public Open Targets Platform
GraphQL API and stores, for every (disease, target) pair, the overall score and the
per-datasource and per-datatype scores.

Sub-commands
------------
list-diseases  Build the disease list (diseases with any evidence of a given datatype,
               default ``genetic_association``) from the official Open Targets bulk
               ``association_by_datatype_direct`` parquet of the SAME release.  The
               GraphQL API has no endpoint that enumerates all diseases by datatype.
fetch          Fetch one chunk of the disease list (SLURM-array friendly) into a
               resumable per-disease cache and write a chunk parquet.
merge          Concatenate chunk parquets into one long table.

Outputs (long format)
---------------------
<out>                       diseaseId, targetId, datasourceId, score
<out stem>_datatype.parquet diseaseId, targetId, datatypeId, score
<out stem>_overall.parquet  diseaseId, targetId, score

Cache layout: <cache-dir>/<release>/<direct|indirect>/<diseaseId>.json.gz
Each cache file is written atomically only after all pages of a disease are retrieved,
so an interrupted job simply resumes with the next uncached disease.

Dependencies: requirements.txt plus ``requests`` (``.venv/bin/pip install requests``), needed only by
``fetch``; ``list-diseases``, ``merge`` and ``--help`` run without it.  Not run by run_pipeline.sh: the
shipped Figure 1 tables (data/README.md, Figure 1) are the record; this script is the optional
first step of rebuilding them (then ``python src/disatlas/fig1_open_targets.py tables``).

Example
-------
python scripts/fetch_open_targets.py fetch --disease-list diseases.txt --out chunk_000.parquet \
    --chunk 0/50 --cache-dir cache
"""
from __future__ import annotations

import argparse
import gzip
import json
import logging
import os
import random
import sys
import time
from pathlib import Path

import pandas as pd

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"
USER_AGENT = "open-targets-association-fetcher/1.0 (+research use)"

META_QUERY = "{ meta { dataVersion { year month iteration } apiVersion { x y z } } }"

ASSOC_QUERY = """
query assoc($efoId: String!, $index: Int!, $size: Int!, $indirect: Boolean!) {
  disease(efoId: $efoId) {
    id
    associatedTargets(page: {index: $index, size: $size}, enableIndirect: $indirect) {
      count
      rows {
        target { id }
        score
        datasourceScores { id score }
        datatypeScores { id score }
      }
    }
  }
}
"""

log = logging.getLogger("fetch_open_targets")

# Data root: --data-root  >  $OT_DATA_ROOT  >  "data_root" in fetch_config.json next to
# this script  >  ./open_targets_data.  Relative --out / --cache-dir / bulk paths are
# resolved against the data root, so a public release only needs to edit the config.
CONFIG_FILE = Path(__file__).with_name("fetch_config.json")


def default_data_root() -> Path:
    if os.environ.get("OT_DATA_ROOT"):
        return Path(os.environ["OT_DATA_ROOT"])
    if CONFIG_FILE.exists():
        root = json.loads(CONFIG_FILE.read_text()).get("data_root")
        if root:
            return Path(root)
    return Path("open_targets_data")


def _under(root: Path, p: str | os.PathLike) -> Path:
    p = Path(p)
    return p if p.is_absolute() else root / p


# --------------------------------------------------------------------------- HTTP
class Client:
    """GraphQL client with polite rate limiting and exponential backoff."""

    def __init__(self, url: str = API_URL, min_interval: float = 0.25,
                 max_retries: int = 8, timeout: float = 120.0):
        self.url = url
        self.min_interval = min_interval
        self.max_retries = max_retries
        self.timeout = timeout
        self._last = 0.0
        import requests  # optional dependency, only for fetch
        self._requests = requests
        self.session = requests.Session()
        self.session.headers.update({"Content-Type": "application/json",
                                     "User-Agent": USER_AGENT})
        self.n_requests = 0
        self.n_retries = 0

    def _wait(self) -> None:
        dt = time.monotonic() - self._last
        if dt < self.min_interval:
            time.sleep(self.min_interval - dt)
        self._last = time.monotonic()

    def query(self, query: str, variables: dict | None = None) -> dict:
        payload = {"query": query, "variables": variables or {}}
        for attempt in range(self.max_retries + 1):
            self._wait()
            self.n_requests += 1
            try:
                r = self.session.post(self.url, json=payload, timeout=self.timeout)
                if r.status_code == 200:
                    body = r.json()
                    if body.get("errors") and not body.get("data"):
                        raise RuntimeError(f"GraphQL error: {body['errors']}")
                    return body["data"]
                retryable = r.status_code in (408, 429, 500, 502, 503, 504)
                if not retryable:
                    r.raise_for_status()
                retry_after = r.headers.get("Retry-After")
                err = f"HTTP {r.status_code}"
            except (self._requests.ConnectionError, self._requests.Timeout, ValueError) as e:
                retry_after, err = None, repr(e)
            if attempt == self.max_retries:
                raise RuntimeError(f"giving up after {attempt + 1} attempts: {err}")
            delay = float(retry_after) if retry_after and retry_after.isdigit() \
                else min(300.0, 2.0 ** attempt) * (1 + random.random() * 0.25)
            self.n_retries += 1
            log.warning("request failed (%s); retry %d in %.1fs", err, attempt + 1, delay)
            time.sleep(delay)
        raise AssertionError("unreachable")

    def release(self) -> str:
        dv = self.query(META_QUERY)["meta"]["dataVersion"]
        rel = f"{int(dv['year']):02d}.{int(dv['month']):02d}"
        if dv.get("iteration") not in (None, "", "0", 0):
            rel += f".{dv['iteration']}"
        return rel


# --------------------------------------------------------------------------- fetch
def fetch_disease(client: Client, efo_id: str, page_size: int, indirect: bool) -> dict:
    rows, index, count = [], 0, None
    while True:
        data = client.query(ASSOC_QUERY, {"efoId": efo_id, "index": index,
                                          "size": page_size, "indirect": indirect})
        dis = data.get("disease")
        if dis is None:  # unknown / obsoleted id in this release
            return {"diseaseId": efo_id, "found": False, "count": 0, "rows": []}
        at = dis["associatedTargets"]
        count = at["count"]
        rows.extend(at["rows"])
        index += 1
        if not at["rows"] or len(rows) >= count:
            break
    if len(rows) != count:
        raise RuntimeError(f"{efo_id}: retrieved {len(rows)} rows, API count {count}")
    return {"diseaseId": efo_id, "found": True, "count": count, "rows": rows}


def _cache_path(cache_dir: Path, release: str, indirect: bool, efo_id: str) -> Path:
    return cache_dir / release / ("indirect" if indirect else "direct") / f"{efo_id}.json.gz"


def _write_atomic(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + f".tmp{os.getpid()}")
    with gzip.open(tmp, "wt") as fh:
        json.dump(obj, fh, separators=(",", ":"))
    os.replace(tmp, path)


def records_to_frames(records: list[dict]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ds, dt, ov = [], [], []
    for rec in records:
        d = rec["diseaseId"]
        for row in rec["rows"]:
            t = row["target"]["id"]
            ov.append((d, t, row["score"]))
            ds.extend((d, t, s["id"], s["score"]) for s in row["datasourceScores"])
            dt.extend((d, t, s["id"], s["score"]) for s in row["datatypeScores"])
    return (pd.DataFrame(ds, columns=["diseaseId", "targetId", "datasourceId", "score"]),
            pd.DataFrame(dt, columns=["diseaseId", "targetId", "datatypeId", "score"]),
            pd.DataFrame(ov, columns=["diseaseId", "targetId", "score"]))


def _sibling(out: Path, tag: str) -> Path:
    return out.with_name(f"{out.stem}_{tag}{out.suffix}")


def read_disease_list(path: Path) -> list[str]:
    if path.suffix == ".parquet":
        ids = pd.read_parquet(path).iloc[:, 0].astype(str).tolist()
    else:
        ids = [ln.strip().split(",")[0] for ln in path.read_text().splitlines()]
        ids = [i for i in ids if i and i not in ("diseaseId", "id")]
    return sorted(set(ids))


def parse_chunk(spec: str | None) -> tuple[int, int]:
    if spec is None:
        if "SLURM_ARRAY_TASK_ID" in os.environ and "SLURM_ARRAY_TASK_COUNT" in os.environ:
            return int(os.environ["SLURM_ARRAY_TASK_ID"]), int(os.environ["SLURM_ARRAY_TASK_COUNT"])
        return 0, 1
    i, n = (int(x) for x in spec.split("/"))
    if not 0 <= i < n:
        raise ValueError(f"bad --chunk {spec}")
    return i, n


def cmd_fetch(a: argparse.Namespace) -> None:
    client = Client(min_interval=a.min_interval, max_retries=a.max_retries)
    release = a.release or client.release()
    root = Path(a.data_root)
    ids = read_disease_list(_under(root, a.disease_list))
    i, n = parse_chunk(a.chunk)
    mine = ids[i::n]
    if a.limit:
        mine = mine[: a.limit]
    log.info("release %s | chunk %d/%d | %d of %d diseases", release, i, n, len(mine), len(ids))
    cache_dir = _under(root, a.cache_dir)
    records, t0 = [], time.time()
    for k, efo in enumerate(mine, 1):
        p = _cache_path(cache_dir, release, a.indirect, efo)
        if p.exists():
            with gzip.open(p, "rt") as fh:
                rec = json.load(fh)
        else:
            rec = fetch_disease(client, efo, a.page_size, a.indirect)
            rec.update(release=release, indirect=a.indirect)
            _write_atomic(p, rec)
        records.append(rec)
        if k % 50 == 0 or k == len(mine):
            log.info("%d/%d diseases | %d requests, %d retries | %.0fs",
                     k, len(mine), client.n_requests, client.n_retries, time.time() - t0)
    ds, dt, ov = records_to_frames(records)
    out = _under(root, a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    ds.to_parquet(out, index=False)
    dt.to_parquet(_sibling(out, "datatype"), index=False)
    ov.to_parquet(_sibling(out, "overall"), index=False)
    missing = [r["diseaseId"] for r in records if not r["found"]]
    json.dump({"release": release, "chunk": f"{i}/{n}", "n_diseases": len(mine),
               "n_not_found": len(missing), "not_found": missing,
               "n_pairs": len(ov), "n_datasource_rows": len(ds),
               "requests": client.n_requests, "retries": client.n_retries,
               "seconds": round(time.time() - t0, 1)},
              open(_sibling(out, "summary").with_suffix(".json"), "w"), indent=1)


def cmd_list_diseases(a: argparse.Namespace) -> None:
    """Diseases with >0 direct association score for any of the requested datatypes.

    Supports both bulk schemas: (diseaseId, datatypeId, score) [<=24.x] and
    (diseaseId, aggregationType, aggregationValue, associationScore) [>=25.x].
    """
    import pyarrow.parquet as pq
    root = Path(a.data_root)
    a.bulk_datatype_dir = _under(root, a.bulk_datatype_dir)
    names = set(pq.read_schema(next(Path(a.bulk_datatype_dir).rglob("*.parquet"))).names)
    if "aggregationValue" in names:
        df = pd.read_parquet(a.bulk_datatype_dir,
                             columns=["diseaseId", "aggregationValue", "associationScore"])
        df.columns = ["diseaseId", "datatypeId", "score"]
    else:
        df = pd.read_parquet(a.bulk_datatype_dir, columns=["diseaseId", "datatypeId", "score"])
    datatypes = [d.strip() for d in a.datatype.split(",")]
    sel = df[df["datatypeId"].isin(datatypes) & (df["score"] > 0)]
    for d in datatypes:
        log.info("%s: %d diseases", d, sel.loc[sel.datatypeId == d, "diseaseId"].nunique())
    ids = sorted(sel["diseaseId"].unique())
    _under(root, a.out).write_text("\n".join(ids) + "\n")
    log.info("%d diseases with evidence in any of %s", len(ids), datatypes)


def cmd_merge(a: argparse.Namespace) -> None:
    root = Path(a.data_root)
    out = _under(root, a.out)
    for tag in (None, "datatype", "overall"):
        parts = sorted(p for p in _under(root, a.chunk_dir).glob(a.pattern)
                       if (tag is None and not p.stem.endswith(("_datatype", "_overall")))
                       or (tag is not None and p.stem.endswith(f"_{tag}")))
        df = pd.concat([pd.read_parquet(p) for p in parts], ignore_index=True)
        target = out if tag is None else _sibling(out, tag)
        df.to_parquet(target, index=False)
        log.info("%s: %d parts, %d rows, %d diseases", target.name, len(parts), len(df),
                 df["diseaseId"].nunique())


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="fetch one chunk of diseases from the GraphQL API")
    f.add_argument("--disease-list", required=True, help="txt (one id per line) or parquet")
    f.add_argument("--out", required=True, help="output parquet (datasource-level long table)")
    f.add_argument("--chunk", default=None, help="'i/N' (default: SLURM array env, else 0/1)")
    f.add_argument("--cache-dir", default="ot_cache")
    f.add_argument("--release", default=None, help="override; default: queried from API meta")
    f.add_argument("--page-size", type=int, default=3000)
    f.add_argument("--indirect", action="store_true", help="include indirect (ontology-propagated) evidence")
    f.add_argument("--min-interval", type=float, default=0.25, help="seconds between requests")
    f.add_argument("--max-retries", type=int, default=8)
    f.add_argument("--limit", type=int, default=0, help="debug: only first N diseases of chunk")
    f.set_defaults(func=cmd_fetch)

    l = sub.add_parser("list-diseases", help="disease list from bulk association_by_datatype_direct")
    l.add_argument("--bulk-datatype-dir", default="association_by_datatype_direct",
                   help="bulk parquet dir (relative to --data-root)")
    l.add_argument("--datatype", default="genetic_association",
                   help="comma-separated datatype ids, e.g. genetic_association,genetic_literature")
    l.add_argument("--out", required=True)
    l.set_defaults(func=cmd_list_diseases)

    m = sub.add_parser("merge", help="merge chunk parquets")
    m.add_argument("--chunk-dir", required=True)
    m.add_argument("--pattern", default="chunk_*.parquet")
    m.add_argument("--out", required=True)
    m.set_defaults(func=cmd_merge)

    for p in (f, l, m):
        p.add_argument("--data-root", default=str(default_data_root()),
                       help="root for relative paths (default: $OT_DATA_ROOT, fetch_config.json, "
                            "or ./open_targets_data)")
    a = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        stream=sys.stderr)
    a.func(a)


if __name__ == "__main__":
    main()
