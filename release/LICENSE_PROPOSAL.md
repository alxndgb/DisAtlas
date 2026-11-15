# Code licence proposal

**Decision (author, 2026-11-15): MIT.** `LICENSE` is committed with the standard MIT text; copyright holder: The DisAtlas
authors (author decision). The proposal below is kept as the record of the options considered.

Original proposal (before the decision):

1. **Code (src/, scripts/, figures/, tables/, tests/, notebooks/, run_pipeline.sh): Apache License 2.0.** Permissive like MIT,
   with an explicit patent grant and a NOTICE mechanism, which corporate open-source policies commonly prefer. MIT is the
   alternative if a shorter licence is wanted. Either requires sign-off from the employer's open-source review.
2. **Data and model (data/, Zenodo record): CC BY 4.0** for author-generated content, with the third-party and share-alike
   exceptions (ChEMBL-derived files keep CC BY-SA 3.0; HPA is CC BY 4.0).
3. Once chosen: add `LICENSE` with the full licence text, a `NOTICE` listing the third-party data attributions from the authors' licence table
   (needed for CC BY sources: Orphanet, GO, Mondo, HPO), set the licence in the package metadata if one is added (the repository has no pyproject.toml), set the Zenodo record licence, and
   remove the TBD marker in README.md.
4. **Unconfirmed data licences** (authors' internal licence table, not in the repository): every data row is still marked
   for author confirmation. Set by the author on 2026-11-15 (not verified on the sources): HPO licence for HPO; CC BY-SA 3.0
   for the Zenodo `event_time/` inputs (ChEMBL share-alike); CC BY 4.0 for the derived salted edges of the training graph,
   the BIKG HABCD subset and the restricted-source tokens, with data-owner clearance still pending.
