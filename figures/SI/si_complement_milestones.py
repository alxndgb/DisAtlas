"""SI figure: clinical milestones of complement inhibitors (C5 and other complement targets), per drug x indication.
Inputs (mechanism_timelines/tables): mechanism_timeline_approved_indications.tsv (approved rows: first trial, first Phase 3,
first trial publication, FDA approval), mechanism_timeline_failed_ongoing_indications.tsv (not approved: failed = Phase 2+
stopped for efficacy; ongoing kept only with an active Phase 3, the source figure's rule), fda_supplements_raw.tsv
(eculizumab PNH original approval, BLA125166 ORIG-1, which the indication table records only as pre-SPL).
Avacopan trial dates come from ClinicalTrials.gov (SI_complement_trial_sources.tsv), because the matched trials post-date its approval:
first registered patient trial NCT01363388 (Phase 2, ANCA-associated vasculitis), Phase 3 NCT02994927 (ADVOCATE), for both GPA and MPA.
When the first registered trial is itself the Phase 3, a ring around the triangle marks it.
Any remaining trial date that post-dates the FDA approval is removed (unsourced).
Every plotted date is written with its source to SI_complement_milestones_dates.tsv.
Source checks (si_complement_milestones_checks: openFDA label + Drugs@FDA, ClinicalTrials.gov, PubMed, fetched records kept in
fda_ctgov_pubmed_records.json): narsoplimab's approved indication is TA-TMA on its label (BLA761152), not aHUS; its row is relabelled
with the TA-TMA pivotal study (NCT02222545, Phase 2) and its first publication; the aHUS Phase 3 (sponsor-terminated) is not
used. A failed row is kept only if ClinicalTrials.gov gives lack of efficacy as the stop reason (SI_complement_failed_stop_reasons.tsv).
Usage: si_complement_milestones.py <src_dir> <out_dir> <stop_reasons.tsv>"""
from pathlib import Path as _Path
import sys as _sys
HERE = _Path(__file__).resolve().parent; REPO = HERE.parents[1]; _sys.path.insert(0, str(REPO / "src"))
from disatlas.sibio_common import DERIVED, ZENODO_SIBIO, RAW  # noqa: E402
import sys, json, numpy as np, pandas as pd, matplotlib as mpl
from matplotlib.lines import Line2D; from matplotlib.patches import Patch, Rectangle
from disatlas import figstyle_si as fs
SRC, OUT = DERIVED / "si_mechanism_timelines", HERE; STOPF = SRC / "SI_complement_failed_stop_reasons.tsv"; fs.apply_figure_style()
D = ["first_trial_date", "first_ph3_date", "pubmed_first_date", "fda_date", "first_label_date"]
TA = pd.read_csv(f"{SRC}/mechanism_timeline_approved_indications.tsv", sep="\t", parse_dates=D)
TF = pd.read_csv(f"{SRC}/mechanism_timeline_failed_ongoing_indications.tsv", sep="\t", parse_dates=D[:2])
FS = pd.read_csv(f"{SRC}/fda_supplements_raw.tsv", sep="\t")
cm_ = lambda d: d[d.mechanism_label.str.startswith("Complement")].copy()
A, F = cm_(TA), cm_(TF); F = F[(F.status_class == "failed") | ((F.status_class == "ongoing") & F.first_ph3_date.notna())]
SR = pd.read_csv(STOPF, sep="\t"); keep = set(map(tuple, SR.loc[SR.stopped_for_lack_of_efficacy, ["drug", "indication"]].values))
dropF = F[(F.status_class == "failed") & ~F.apply(lambda r: (r.drug, r.disease_label) in keep, axis=1)]
F = F.drop(dropF.index); assert set(map(tuple, SR[["drug", "indication"]].values)) >= set(map(tuple, dropF[["drug", "disease_label"]].values))
na = A.drug == "narsoplimab"; assert na.sum() == 1 and (A.loc[na, "disease_label"] == "atypical hemolytic-uremic syndrome").all() and (A.loc[na, "fda_sub"] == "BLA761152 ORIG-1").all()
A.loc[na, ["disease_label", "first_trial_nct", "first_trial_date", "first_ph3_nct", "first_ph3_date", "pubmed_first_date"]] = [
    "transplant-associated thrombotic microangiopathy", "NCT02222545", pd.Timestamp("2014-08-18"), np.nan, pd.NaT, pd.Timestamp("2022-08-01")]
o = FS[(FS.application == "BLA125166") & (FS.submission_type == "ORIG") & (FS.submission_number == 1)]; assert len(o) == 1
ecu = (A.drug == "eculizumab") & (A.disease_label == "paroxysmal nocturnal hemoglobinuria"); assert A.loc[ecu, "fda_date"].isna().all()
A.loc[ecu, "fda_date"] = pd.to_datetime(str(o.approval_date.iloc[0])); A.loc[ecu, "fda_src"] = "openFDA BLA125166 ORIG-1"
A["fda_src"] = A.fda_src.fillna("openFDA " + A.fda_sub.astype(str))
av = A.drug == "avacopan"; assert set(A.loc[av, "disease_label"]) == {"granulomatosis with polyangiitis", "microscopic polyangiitis"}
A.loc[av, ["first_trial_nct", "first_trial_date", "first_ph3_nct", "first_ph3_date"]] = ["NCT01363388", pd.Timestamp("2011-05-26"), "NCT02994927", pd.Timestamp("2016-12-11")]
bad = A.first_trial_date > A.fda_date; A.loc[bad, ["first_trial_date", "first_ph3_date", "first_trial_nct", "first_ph3_nct"]] = np.nan
yf = lambda d: np.nan if pd.isna(d) else d.year + (d.dayofyear - .5) / 365.25
order = A.groupby("drug").fda_date.min().sort_values().index
rows, trace = [], []
for dr in order:
    a = A[A.drug == dr].assign(k=lambda x: x.role.map({"first_approved": 0, "expansion": 1})).sort_values(["k", "fda_date"])
    for r in a.itertuples(): rows.append(dict(drug=dr, dis=r.disease_label, cls="approved", ft=r.first_trial_date, p3=r.first_ph3_date, pub=r.pubmed_first_date, fda=r.fda_date))
    for r in a.itertuples():
        for ev, d, s in [("first trial", r.first_trial_date, f"ClinicalTrials.gov {r.first_trial_nct}"), ("first Phase 3", r.first_ph3_date, f"ClinicalTrials.gov {r.first_ph3_nct}"),
                         ("first trial publication", r.pubmed_first_date, "PubMed"), ("FDA approval", r.fda_date, r.fda_src)]:
            if pd.notna(d): trace.append(dict(drug=dr, indication=r.disease_label, status="approved", event=ev, date=d.date(), source=s))
    for r in F[F.drug == dr].sort_values("first_trial_date").itertuples():
        rows.append(dict(drug=dr, dis=r.disease_label, cls=r.status_class, ft=r.first_trial_date, p3=r.first_ph3_date, pub=pd.NaT, fda=pd.NaT))
        for ev, d, s in [("first trial", r.first_trial_date, f"ClinicalTrials.gov {r.first_trial_nct}"), ("first Phase 3", r.first_ph3_date, f"ClinicalTrials.gov {r.first_ph3_nct}")]:
            if pd.notna(d): trace.append(dict(drug=dr, indication=r.disease_label, status=r.status_class, event=ev, date=d.date(), source=s))
R = pd.DataFrame(rows); n = len(R); HL = {("eculizumab", "myasthenia gravis"), ("eculizumab", "neuromyelitis optica")}
assert all(((R.drug == d) & (R.dis == x)).sum() == 1 for d, x in HL)
BLUE, GREY, END = "#2166ac", "0.6", 2026.9
fig, ax = mpl.pyplot.subplots(figsize=(7.0, 0.165 * n + 1.6))
for i, r in R.iterrows():
    hl = (r.drug, r.dis) in HL
    if hl: ax.add_patch(Rectangle((1994.6, i - .5), END - 1994.6, 1, fc="#fde9c8", ec="none", zorder=0))
    col = BLUE if r.cls == "approved" else GREY; s0 = yf(r.ft)
    if r.cls == "approved":
        if pd.notna(s0): ax.plot([s0, yf(r.fda)], [i, i], color=col, lw=1.6, zorder=2)
        ax.plot(yf(r.fda), i, "*", ms=7, color=col, mec="k", mew=.4, zorder=5)
        if pd.notna(r.pub): ax.plot(yf(r.pub), i, "D", ms=3.2, mfc="white", mec="k", mew=.6, zorder=6)
    else:
        ax.plot([s0, END], [i, i], color=col, lw=1.6, ls="-" if r.cls == "failed" else (0, (2, 1.5)), zorder=2)
    same = pd.notna(s0) and pd.notna(r.p3) and abs(s0 - yf(r.p3)) < 1e-9
    if same: ax.plot(s0, i, "o", ms=7.2, mfc="white", mec=col, mew=.9, zorder=3)
    elif pd.notna(s0): ax.plot(s0, i, "o", ms=3.8, mfc="white", mec=col, mew=1, zorder=4)
    if pd.notna(r.p3): ax.plot(yf(r.p3), i, "^", ms=4, color=col, zorder=4)
ax.set_yticks(range(n)); ax.set_yticklabels([f"{r.drug} – {r.dis}" for r in R.itertuples()], fontsize=6)
for t, r in zip(ax.get_yticklabels(), R.itertuples()):
    if r.cls != "approved": t.set_color("0.4")
    if (r.drug, r.dis) in HL: t.set_fontweight("bold")
ax.set_ylim(n - .4, -.6); ax.set_xlim(1994.6, END); ax.set_xticks(range(1995, 2030, 5)); ax.set_xlabel("Year")
ax.spines[["top", "right", "left"]].set_visible(False); ax.tick_params(axis="y", length=0)
ax.set_title("Clinical milestones of complement inhibitors by drug and indication", fontsize=8, loc="left")
h = [Line2D([], [], marker="o", ls="", mfc="white", mec="k", ms=4, label="first registered trial"),
     Line2D([], [], marker="^", ls="", color="k", ms=4, label="first Phase 3 trial"),
     (Line2D([], [], marker="o", ls="", mfc="white", mec="k", ms=7.2, mew=.9), Line2D([], [], marker="^", ls="", color="k", ms=4)),
     Line2D([], [], marker="D", ls="", mfc="white", mec="k", ms=3.2, label="first trial publication"),
     Line2D([], [], marker="*", ls="", color=BLUE, mec="k", ms=7, label="FDA approval"),
     Line2D([], [], color=BLUE, lw=1.6, label="approved indication"), Line2D([], [], color=GREY, lw=1.6, label="not approved: stopped for lack of efficacy"),
     Line2D([], [], color=GREY, lw=1.6, ls=(0, (2, 1.5)), label="not approved: Phase 3 ongoing"), Patch(fc="#fde9c8", label="eculizumab case discussed in the text")]
L = ["first registered trial", "first Phase 3 trial", "first registered trial is the Phase 3", "first trial publication", "FDA approval", "approved indication",
     "not approved: stopped for lack of efficacy", "not approved: Phase 3 ongoing", "eculizumab case discussed in the text"]
fig.legend(handles=h, labels=L, loc="lower center", ncol=3, fontsize=6, frameon=False, bbox_to_anchor=(.5, -.01))
fig.subplots_adjust(bottom=1.0 / (0.165 * n + 1.6), left=.36, right=.98, top=1 - .35 / (0.165 * n + 1.6))
pd.DataFrame(trace).to_csv(f"{OUT}/SI_complement_milestones_dates.tsv", sep="\t", index=False)
N = dict(rows=n, approved=int((R.cls == "approved").sum()), failed=int((R.cls == "failed").sum()), ongoing=int((R.cls == "ongoing").sum()),
         drugs=len(order), failed_rows_dropped_stop_reason_not_efficacy=dropF[["drug", "disease_label"]].values.tolist(), trial_markers_omitted=A.loc[bad, ["drug", "disease_label"]].values.tolist(), dated_events=len(trace))
N["png_bytes"] = fs.save(fig, OUT, "SI_complement_milestones"); N["overlaps"] = fs.overlap_check(fig)
json.dump(N, open(f"{OUT}/SI_complement_milestones_numbers.json", "w"), indent=1, default=str); print(json.dumps(N, default=str))
