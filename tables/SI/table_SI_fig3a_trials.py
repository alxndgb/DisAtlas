"""Appendix, "Trials behind the Figure 3a example" -- Table tbl:SI_fig3a_trials.

Lists every trial behind the ten disease-drug lines of Figure 3a.
  tables/inputs/si_fig3a_trials.tsv     one row per (trial, disease, drug): trial_key as in
                                        data/00_raw/trials/trial_triples.parquet, public registry ID
                                        (empty when the trial has none), phase as printed,
                                        ClinicalTrials.gov phase, title
  data/00_raw/trials/trial_triples.parquet   drugs in trial, start year, highest phase
  data/00_raw/diseases, data/00_raw/drugs     disease and drug labels
"""
import sys
from pathlib import Path

import pandas as pd

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tablib import INPUTS, RAW, write_and_check  # noqa: E402

PHASE_ORDER = {"IV": 4, "III": 3, "II/III": 2.5, "II": 2, "I/II": 1.5, "I": 1, "OTHER": 0}
PHASE_INT = {"IV": 4, "III": 3, "II": 2, "I/II": 2, "I": 1}

t = pd.read_csv(INPUTS / "si_fig3a_trials.tsv", sep="\t", dtype=str, keep_default_na=False, na_values=[""])
tt = pd.read_parquet(RAW / "trials" / "trial_triples.parquet")
tt = tt[tt.disease_id.isin(t.disease_id) & tt.drug_id.isin(t.drug_id)]
t = t.merge(tt[["trial_key", "disease_id", "drug_id", "n_drugs_in_trial", "max_phase", "year"]],
            on=["trial_key", "disease_id", "drug_id"], how="left", validate="1:1")
assert t.n_drugs_in_trial.notna().all(), "trial not in trial_triples"
n_line = tt.groupby(["disease_id", "drug_id"]).trial_key.nunique()
assert n_line.to_dict() == t.groupby(["disease_id", "drug_id"]).trial_key.nunique().to_dict()
ph = t.phase.map(PHASE_INT)
assert (ph[ph.notna()] == t.max_phase[ph.notna()]).all(), "phase differs from trial_triples"

dis = pd.read_parquet(RAW / "diseases" / "diseases.parquet").set_index("disease_id").label
drg = pd.read_parquet(RAW / "drugs" / "drugs.parquet").set_index("drug_id").name
t["disease"] = t.disease_id.map(dis)
t["drug"] = t.drug_id.map(drg).str.lower()
t["_o"] = t.phase.map(PHASE_ORDER)
t = t.sort_values(["disease", "drug", "_o", "year", "public_id"], ascending=[True, True, False, True, True],
                  na_position="last", kind="stable").reset_index(drop=True)
differ = t.ctgov_phase.notna() & (t.ctgov_phase != t.phase)
cell_phase = ["other" if p == "OTHER" else p for p in t.phase]
df = pd.DataFrame({
    "Disease": t.disease, "Drug": t.drug, "Trial": t.public_id.fillna("–"),
    "Phase": [p + "†" if d else p for p, d in zip(cell_phase, differ)],
    "Start": ["–" if pd.isna(y) else str(int(y)) for y in t.year],
    "Drugs in trial": t.n_drugs_in_trial.astype(int).astype(str), "Title": t.title})

n, n_id, n_mono = len(t), int(t.public_id.notna().sum()), int((t.n_drugs_in_trial == 1).sum())
n_diff = int(differ.sum())
n_na = int((t.ctgov_phase[differ] == "NA").sum())
assert set(t.ctgov_phase[differ]) <= {"NA", "I"}
assert (n, n_id, n_mono, n_diff, n_na, n_diff - n_na) == (90, 72, 61, 7, 6, 1)

EXPECTED = [
    ['Alzheimer disease', 'creatine', 'NCT05383833', 'IV†', '2022', '1'],
    ['Alzheimer disease', 'pimavanserin', 'NCT03325556', 'III', '2017', '1'],
    ['Alzheimer disease', 'pimavanserin', 'NCT03575052', 'III', '2018', '1'],
    ['Alzheimer disease', 'pimavanserin', 'NCT03623321', 'III', '2018', '1'],
    ['Alzheimer disease', 'pimavanserin', '–', 'II', '2010', '1'],
    ['Alzheimer disease', 'pimavanserin', 'NCT02035553', 'II', '2013', '1'],
    ['Alzheimer disease', 'pimavanserin', 'NCT02992132', 'II', '2016', '1'],
    ['Alzheimer disease', 'pimavanserin', 'NCT03118947', 'II', '2017', '1'],
    ['Alzheimer disease', 'prazosin', 'NCT01126099', 'IV†', '2010', '1'],
    ['Alzheimer disease', 'prazosin', 'NCT00161473', 'II†', '2001', '1'],
    ['Alzheimer disease', 'prazosin', 'NCT03710642', 'II', '2018', '1'],
    ['Alzheimer disease', 'prazosin', '–', 'II', '–', '1'],
    ['Alzheimer disease', 'trazodone', 'NCT05209035', 'IV', '2022', '1'],
    ['Alzheimer disease', 'trazodone', 'NCT00898807', 'III', '2009', '6'],
    ['Alzheimer disease', 'trazodone', 'NCT01142258', 'III', '2010', '1'],
    ['Alzheimer disease', 'trazodone', 'NCT00000179', 'III', '–', '3'],
    ['Alzheimer disease', 'trazodone', 'NCT01584440', 'II', '2012', '5'],
    ['Alzheimer disease', 'trazodone', 'NCT05282550', 'II', '2023', '1'],
    ['Alzheimer disease', 'trazodone', '–', 'other', '–', '11'],
    ['Lewy body dementia', 'pimavanserin', 'NCT05590637', 'IV', '2022', '2'],
    ['Lewy body dementia', 'pimavanserin', 'NCT05357612', 'IV', '2023', '1'],
    ['Lewy body dementia', 'pimavanserin', 'NCT03325556', 'III', '2017', '1'],
    ['Lewy body dementia', 'pimavanserin', 'NCT03575052', 'III', '2018', '1'],
    ['Lewy body dementia', 'pimavanserin', 'NCT03623321', 'III', '2018', '1'],
    ['Parkinson disease', 'creatine', '–', 'IV', '–', '2'],
    ['Parkinson disease', 'creatine', '–', 'II', '–', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT03152292', 'IV', '2017', '2'],
    ['Parkinson disease', 'pimavanserin', 'NCT04292223', 'IV', '2020', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT04373317', 'IV', '2022', '2'],
    ['Parkinson disease', 'pimavanserin', 'NCT05590637', 'IV', '2022', '2'],
    ['Parkinson disease', 'pimavanserin', 'NCT05357612', 'IV', '2023', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT00477672', 'III', '2007', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT00550238', 'III', '2007', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT00658567', 'III', '2008', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT01174004', 'III', '2010', '1'],
    ['Parkinson disease', 'pimavanserin', '–', 'III', '2013', '1'],
    ['Parkinson disease', 'pimavanserin', 'EudraCT 2016-003512-12', 'III', '2017', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT03325556', 'III', '2017', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT03575052', 'III', '2018', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT03623321', 'III', '2018', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT06068465', 'III', '2023', '2'],
    ['Parkinson disease', 'pimavanserin', 'CTR20251027', 'III', '2025', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20252689', 'III', '2025', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT02762591', 'III', '–', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT00086294', 'II', '2004', '2'],
    ['Parkinson disease', 'pimavanserin', 'NCT00087542', 'II', '2004', '2'],
    ['Parkinson disease', 'pimavanserin', 'NCT01518309', 'II', '2004', '2'],
    ['Parkinson disease', 'pimavanserin', 'NCT03482882', 'II', '2018', '4'],
    ['Parkinson disease', 'pimavanserin', 'NCT04164758', 'II', '2019', '3'],
    ['Parkinson disease', 'pimavanserin', 'NCT03947216', 'II', '2020', '1'],
    ['Parkinson disease', 'pimavanserin', '–', 'II', '–', '1'],
    ['Parkinson disease', 'pimavanserin', '–', 'I/II', '–', '1'],
    ['Parkinson disease', 'pimavanserin', '–', 'I', '2003', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20220014', 'I', '2022', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20222082', 'I', '2022', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20222456', 'I', '2022', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20222909', 'I', '2023', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20230288', 'I', '2023', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20231339', 'I', '2023', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20232504', 'I', '2023', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT05796167', 'I', '2023', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20240780', 'I', '2024', '1'],
    ['Parkinson disease', 'pimavanserin', 'CTR20241274', 'I', '2024', '1'],
    ['Parkinson disease', 'pimavanserin', 'NCT06450184', 'I', '2024', '1'],
    ['Parkinson disease', 'pimavanserin', '–', 'I', '–', '1'],
    ['Parkinson disease', 'trazodone', 'IRCT20170821035819N2', 'III', '2020', '3'],
    ['Parkinson disease', 'trazodone', '–', 'II', '2007', '6'],
    ['Parkinson disease', 'trazodone', 'RBR-4b32m3', 'II', '2020', '1'],
    ['Parkinson disease', 'trazodone', '–', 'II', '–', '4'],
    ['alcohol abuse', 'prazosin', 'NCT00167687', 'IV', '2005', '2'],
    ['alcohol abuse', 'prazosin', 'NCT02226367', 'IV†', '2015', '2'],
    ['alcohol abuse', 'prazosin', '–', 'IV', '–', '1'],
    ['alcohol abuse', 'prazosin', 'NCT00744055', 'III', '2009', '1'],
    ['alcohol abuse', 'prazosin', 'NCT00240227', 'II†', '2004', '1'],
    ['alcohol abuse', 'prazosin', 'NCT00762710', 'II', '2008', '1'],
    ['alcohol abuse', 'prazosin', 'NCT01518972', 'II', '2009', '1'],
    ['alcohol abuse', 'prazosin', 'NCT01916941', 'II†', '2013', '1'],
    ['alcohol abuse', 'prazosin', 'NCT02193256', 'II†', '2014', '2'],
    ['alcohol abuse', 'prazosin', 'NCT02322047', 'II', '2015', '2'],
    ['alcohol abuse', 'prazosin', 'NCT04793685', 'II', '2021', '1'],
    ['alcohol abuse', 'prazosin', 'NCT00585780', 'I/II', '2009', '2'],
    ['alcohol abuse', 'prazosin', 'NCT02966340', 'I/II', '2016', '1'],
    ['alcohol abuse', 'prazosin', 'NCT06147622', 'I', '2024', '3'],
    ['alcohol abuse', 'trazodone', 'NCT01662297', 'IV', '2012', '2'],
    ['alcohol abuse', 'trazodone', 'CTRI/2025/08/093625', 'IV', '2025', '2'],
    ['alcohol abuse', 'trazodone', '–', 'IV', '–', '1'],
    ['alcohol abuse', 'trazodone', '–', 'IV', '–', '1'],
    ['alcohol abuse', 'trazodone', '–', 'II', '–', '2'],
    ['alcohol abuse', 'trazodone', '–', 'II', '–', '3'],
    ['alcohol abuse', 'trazodone', '–', 'II', '–', '2'],
]
assert df.iloc[:, :6].values.tolist() == [list(r) for r in EXPECTED], "rows differ from the manuscript"

if __name__ == "__main__":
    write_and_check(df, __file__, "tbl_SI_fig3a_trials", None,
                    f"Trials behind the disease–drug lines of Figure 3a. The {n} trials linking the four "
                    f"diseases to the four drugs, as recorded in the trial data used to build IxIDN; {n_mono} "
                    f"tested a single drug. Trial gives the public registry identifier (ClinicalTrials.gov, "
                    f"EudraCT, China Drug Trials, IRCT, ReBEC or CTRI; {n_id} of {n}), or – when the trial has "
                    f"none. Phase is the highest phase in the trial data; † marks the {n_diff} trials whose "
                    f"ClinicalTrials.gov record gives a different phase ({n_na} not applicable, "
                    f"{n_diff - n_na} phase I). Start is the trial start year (– when not recorded).")
    print(f"tbl_SI_fig3a_trials: {n} rows x {df.shape[1]} cols, rows and caption counts equal the manuscript")
