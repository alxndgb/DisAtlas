"""Raw graph layer: released BIKG edge file integrity and the rebuilt Habcd graph signature.

The full rebuild (distances over ~5.5M candidate pairs) is opt-in: DISATLAS_RAWGRAPH=1 (about 1 h, one core)
and needs the term-distance slices in data/zenodo/ontology_term_term_distance.
"""
import hashlib
import json
import os

import pandas as pd
import pytest

from disatlas import paths
from disatlas.graph import raw_build

REL = paths.RAW / "bikg_habcd" / "habcd_bikg_gene_edges.parquet"
PUB = paths.RAW / "public_disease_gene" / "disease_gene.parquet"
SIG = json.load(open(paths.ROOT / "tests" / "raw_graph_signature.json"))


def test_bikg_release_file():
    sha = open(str(REL) + ".sha256").read().split()[0]
    assert hashlib.sha256(REL.read_bytes()).hexdigest() == sha == SIG["bikg_release_sha256"]
    e = pd.read_parquet(REL)
    assert list(e.columns) == ["disease_id", "gene_symbol", "evidence_tier"]
    assert len(e) == 8343 and not e.duplicated(["disease_id", "gene_symbol"]).any()


def test_assembled_edge_tables():
    g = raw_build.gene_edges()
    p = raw_build.phenotype_edges()
    assert len(g) == SIG["n_gene_edges"] and not g.duplicated(["disease_id", "gene_symbol"]).any()
    assert g.evidence_tier.value_counts().sort_index().to_dict() == {int(k): v for k, v in SIG["gene_edges_by_tier"].items()}
    assert len(p) == SIG["n_phenotype_edges"]


@pytest.mark.skipif(os.environ.get("DISATLAS_RAWGRAPH") != "1" or not raw_build.ZEN.is_dir(),
                    reason="set DISATLAS_RAWGRAPH=1 with data/zenodo/ontology_term_term_distance present")
def test_rebuilt_graph_signature():
    e, q, _ = raw_build.build()
    assert q == SIG["quotas"]
    assert (len(e), len(set(e.u) | set(e.v))) == (SIG["n_edges"], SIG["n_nodes"])
    assert raw_build.edge_signature(e) == SIG["edge_sha256"]


def test_public_disease_gene_file():
    e = pd.read_parquet(PUB)
    assert list(e.columns) == ["disease_id", "gene_symbol", "source"] and len(e) == 308
    assert not e.duplicated(["disease_id", "gene_symbol"]).any()
    assert set(e.source) == {"Open Targets"}
    assert not (paths.RAW / "curated_unattributed").exists()


@pytest.mark.skipif(not (paths.ZENODO / "graph" / "graph_edges.parquet").is_file(),
                    reason="data/zenodo/graph not present")
def test_shipped_broad_graph_unchanged():
    e = pd.read_parquet(paths.ZENODO / "graph" / "graph_edges.parquet")
    s = SIG["shipped_broad"]
    assert (len(e), len(set(e.u) | set(e.v))) == (s["n_edges"], s["n_nodes"])
    assert raw_build.edge_signature(e) == s["edge_sha256"]
