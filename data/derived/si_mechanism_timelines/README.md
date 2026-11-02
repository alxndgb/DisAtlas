# Mechanism timelines and complement milestones

Small precomputed tables that cannot be regenerated inside this repository (upstream model runs or benchmark construction). Each line gives the content and the scripts that read the file.

- `SI_complement_failed_stop_reasons.tsv`: stop reasons of failed complement-inhibitor trials. Read by: si_complement_milestones.py.
- `disatlas_slice_scores_5seed.tsv`: five-seed mean matched percentile per pair and slice. Read by: si_mechanism_coverage_slices.py.
- `fda_supplements_raw.tsv`: FDA submission records (public openFDA) used for approval dates. Read by: si_complement_milestones.py.
- `mechanism_class_controls.tsv`: mechanism-class controls: cases and majority-recovered counts per class. Read by: si_bio_gap_numbers.py.
- `mechanism_timeline_approved_indications.tsv`: approved drug-indication rows with milestone dates. Read by: si_complement_milestones.py.
- `mechanism_timeline_failed_ongoing_indications.tsv`: failed and ongoing drug-indication rows. Read by: si_complement_milestones.py.
- `pair_lead_times.tsv`: first crossing slice and decoy crossing rates per pair. Read by: si_mechanism_coverage_slices.py.

Licence: drugs are FDA-approved products and every date or trial identifier comes from public records (openFDA/Drugs@FDA submissions and labels, ClinicalTrials.gov); no licensed-source fields are shipped.
