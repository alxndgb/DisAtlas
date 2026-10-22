"""Step 1 — clinical-trial disease–drug (DD) network and its layer-0/1 hygiene.

Manuscript: Methods, "Indication-based disease network (IxIDN)", layers 0-1
(non-disease indications removed; non-mechanistic drugs removed; trials testing more
than two kept drugs ignored).

Inputs (data/00_raw/): trials/trial_triples.parquet, diseases/diseases.parquet,
drugs/drug_moa_decision.parquet; data/zenodo/00_raw/ontologies/mondo.obo.
Outputs (data/01_rebuilt/):
  dd_nodes.parquet       one row per DD-network node: node_id, kind (disease/drug), label,
                         is_disease_indication, is_rare (diseases), degree, decision (drugs)
  dd_survivors.parquet   one row per trial record surviving layers 0-1, in record order, with
                         passes_size_cap (False = trial tests > 2 kept drugs; used only by
                         the co-investigation veto, which ignores the size cap)
"""
import sys
import os

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ixidn_lib as L  # noqa: E402


def main():
    k = L.KNOBS
    t = L.load_triples()
    dis = L.load_diseases()
    dec = L.load_drug_decisions()
    rare_mondo = L.mondo_rare()

    # ---- nodes ----
    dd_dis = sorted(set(t.disease_id))
    dd_drug = sorted(set(t.drug_id))
    dmeta = dis.set_index("disease_id")
    ind = dmeta.is_disease_indication.reindex(dd_dis)
    assert ind.notna().all(), "every DD-network disease carries is_disease_indication"
    ddeg = t.groupby("disease_id").drug_id.nunique()
    gdeg = t.groupby("drug_id").disease_id.nunique()
    dn = pd.DataFrame({
        "node_id": dd_dis, "kind": "disease",
        "label": dmeta.label.reindex(dd_dis).values,
        "is_disease_indication": ind.astype(bool).values,
        # rare = the node is a MONDO term in the orphanet_rare | gard_rare subset
        "is_rare": [n in rare_mondo for n in dd_dis],
        "degree": ddeg.reindex(dd_dis).values, "decision": None})
    dm = dec.set_index("drug_id")
    gn = pd.DataFrame({
        "node_id": dd_drug, "kind": "drug", "label": None,
        "is_disease_indication": None, "is_rare": None,
        "degree": gdeg.reindex(dd_drug).values,
        "decision": dm.decision.reindex(dd_drug).values})
    nodes = pd.concat([dn, gn], ignore_index=True)
    nodes.to_parquet(L.out_path("dd_nodes.parquet"), index=False)

    # ---- layer 0/1 survivors ----
    keep_dis = set(dn.node_id[dn.is_disease_indication])
    drop_drug = set(gn.node_id[gn.decision.isin(k["drop_decisions"])])
    s = t[t.disease_id.isin(keep_dis) & ~t.drug_id.isin(drop_drug)].copy()
    n_kept = (s.n_drugs_in_trial.astype("Int64") - s.n_excluded_drugs_in_trial.astype("Int64")).clip(lower=1)
    # unknown trial size: kept (cannot be shown to be a basket trial)
    s["passes_size_cap"] = (n_kept.isna() | (n_kept <= k["max_drugs_per_trial"])).astype(bool)
    s = s.reset_index(drop=True)   # record order kept: step 2 dedups per trial, first wins
    s.to_parquet(L.out_path("dd_survivors.parquet"), index=False)
    print(f"DD network: {len(dd_dis)} diseases ({len(keep_dis)} disease indications), "
          f"{len(dd_drug)} drugs ({len(drop_drug)} non-mechanistic); "
          f"{len(s)} surviving trial records, {int(s.passes_size_cap.sum())} within the size cap")


if __name__ == "__main__":
    main()
