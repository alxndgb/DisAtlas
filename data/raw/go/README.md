# data/raw/go

`go_bp_terms.tsv.gz` (term_id, name, namespace, ancestors, organ_level): Gene Ontology terms annotated in `data/raw/specific_channels/gene2go.parquet`
plus their is_a / part_of ancestors, from go-basic.obo of the authors' working tree (release dated 2026-07-26 in the file name). `ancestors` is '|'-joined;
`organ_level` marks terms at or under GO:0048856 (anatomical structure development) or GO:0003008 (system process).
Written by `scripts/maintainer/stage_cluster_theme_inputs.py`; read by `src/disatlas/rebuild/cluster_themes.py` (tests/test_rebuild_cluster_themes.py).
Licence: GO is distributed under CC BY 4.0 (to be confirmed by author against the release page).

