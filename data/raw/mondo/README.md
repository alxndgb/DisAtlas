# MONDO body-system extract (Figure S16 c)

Staged by `scripts/maintainer/stage_mondo_body_system.py` from MONDO release 2026-02-06 (`mondo.obo`); read by
`src/disatlas/rebuild/cluster_mondo_concordance.py` (test: `tests/test_rebuild_cluster_mondo_concordance.py`).

- `body_system_classes.tsv`: the 19 direct is_a children of MONDO:7770006 "disease by body system or component" (id, name).
- `mondo_body_system_13021.tsv`: for each of the 13,021 clustered diseases, the MONDO term it maps to (itself, or the
  unique MONDO term cross-referencing it; empty if none or ambiguous) and its body-system classes ('|'-joined, via is_a
  ancestors). 12,533 diseases map to MONDO; 10,133 have at least one class.

Source: Mondo Disease Ontology, https://mondo.monarchinitiative.org/ . Licence: to be confirmed by author.
