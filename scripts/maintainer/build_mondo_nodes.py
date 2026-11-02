"""Maintainer-only: derive data/derived/si_compositionality/mondo_nodes.tsv from a MONDO release (mondo.obo).

Not needed to reproduce any result: the table ships with the repository and is what si_compositionality_recompute.py
reads when no ontology file is given.  It keeps only what that analysis uses: for each of the 13,021 embedded nodes its
MONDO term (direct, or the unique non-obsolete MONDO term carrying the node id as an xref), the term name, and the term's
is_a parents that are themselves embedded nodes.  The md5 of the source file is written next to it (mondo_nodes.source.json)
so the recompute records the same input checksum as a run from the ontology file.

usage: python scripts/maintainer/build_mondo_nodes.py <mondo.obo> [--model-dir <released model package>]
"""
import argparse, hashlib, json, sys
from pathlib import Path
REPO = Path(__file__).resolve().parents[2]; sys.path.insert(0, str(REPO / "src"))
from disatlas.mondo_nodes import parse_obo, node_table, TABLE, SOURCE  # noqa: E402
from disatlas.sibio_common import model_dir  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("obo"); ap.add_argument("--model-dir")
a = ap.parse_args()
nodes = open(Path(a.model_dir or model_dir()) / "nodes.txt").read().split(); assert len(nodes) == 13021
T = node_table(parse_obo(a.obo), nodes)
T.to_csv(TABLE, sep="\t", index=False, lineterminator="\n")
json.dump({"mondo.obo": hashlib.md5(open(a.obo, "rb").read()).hexdigest()}, open(SOURCE, "w"), indent=1)
open(SOURCE, "a").write("\n")
print(TABLE, len(T), "nodes;", int(T.mondo_id.notna().sum()), "with a MONDO term")
