"""Main Figure 7 and SI event-time figures: simulated prospective (event-time) validation, conditional-logit design.
Main figure (default): a, post hoc case (eculizumab, PNH -> NMO spectrum disorder; DisAtlas rank of NMO among the
held-out candidates of decision year 2007); b, yearly odds ratio per 1 SD of the DisAtlas score that a drug is next trialed
in a candidate disease (conditional logistic regression, one stratum per drug and decision year, popularity spline,
Efron ties, drug-clustered 95% CIs), unadjusted, adjusted and shuffled-score control; c, pooled 2000-2021 odds ratio per
method, adjusted for popularity.
SI figures (--si): hits10 = share of drugs whose next trial ranks in the top 10% of candidates (pre-specified secondary endpoint),
popularity-matched and no matching; weighted = the same share per drug-year, popularity-weighted with exact popularity
bins, with popularity under the same weights and chance.
Inputs (data/derived/fig7_prospective/, staged by scripts/maintainer/stage_fig7_or.py; the conditional-logit fits were run on
the licensed per-candidate trial rows with scripts/maintainer/fig7_clogit_fit.R and ship as or_clogit.tsv):
  or_clogit.tsv                                  odds ratios (sets 1996-2021 and pooled; models M0, M2, M2_<method>, PERM_M2, H2H)
  nmo_case_clogit.json, nmo_case_scatter.tsv     panel a case (2007 held-out candidates)
  merge_and_matching_census.tsv, ../si_event_time/point_estimates_test.tsv   test-pool census and point estimates
  abs_percentiles.tsv, auc_per_drug.tsv, auc_check_yearly_box.tsv, e2_per_drug_from2000.tsv, weighted_yearly.tsv
  nmosd_timeline.tsv, fig7_nmo_ms_yearly_ens5_test.tsv, nmo_ms_yearly_all_methods.tsv
Literature dates (PNH pilot ~2002, phase 3 2004, AQP4 antibody 2004, NMO criteria 2006, NMOSD criteria 2015) are drawn
as constants and asserted against the sourced milestone rows of nmosd_timeline.tsv.
    python figures/main/fig7_prospective.py [OUTDIR]                 -> fig7_prospective.{png,pdf}, fig7_prospective_values.tsv
    python figures/main/fig7_prospective.py --si hits10 [OUTDIR]     -> SI_event_time_hits10.{png,pdf}
    python figures/main/fig7_prospective.py --si weighted [OUTDIR]   -> SI_event_time_weighted.{png,pdf}
Decision years from 2000 are plotted (earlier years have small comparison sets; they are in the SI tables).
"""
import os, sys
import numpy as np, pandas as pd, matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes (must precede savefig)
SI = None
if "--si" in sys.argv:
    k = sys.argv.index("--si"); SI = sys.argv[k + 1]; del sys.argv[k:k + 2]
assert SI in (None, "hits10", "weighted"), SI
WEIGHTED = SI == "weighted"
PANEL = {None: "or", "hits10": "hits", "weighted": "p1"}[SI]
BONLY = SI is not None                 # SI figures: count strip + panel b only
WV, WU = "wex", "drug_year"            # weighted SI figure: exact popularity bins, pre-specified E2 per drug-year
OUT = sys.argv[1] if len(sys.argv) > 1 else (HERE if SI is None else os.path.join(ROOT, "figures", "SI"))
ET = NMO = FR = OD = os.path.join(ROOT, "data", "derived", "fig7_prospective")
YEARLY = os.path.join(ET, "fig7_nmo_ms_yearly_ens5_test.tsv")
YALL = os.path.join(ET, "nmo_ms_yearly_all_methods.tsv")   # all methods, all pools; baseline rank in panel a
CD='#1F4E99'; PUf='#A893C4'; GR='#555555'; G='#9E9E9E'; ANC='#B35806'
ordn = lambda k: {1: "1st", 2: "2nd", 3: "3rd"}.get(k, f"{k}th")
def frac_year(d):
    d = pd.Timestamp(d); return d.year + (d.dayofyear - 1) / 365.25

# ---- SI panel years and drug counts (test pool census)
CE = pd.read_csv(f"{ET}/merge_and_matching_census.tsv", sep="\t"); CE = CE[CE.pool == "test"]
cen = CE[CE.scope.astype(str).str.isdigit()].copy(); cen["yr"] = cen.scope.astype(int); cen = cen.sort_values("yr")
MISSING = sorted(set(range(cen.yr.min(), cen.yr.max() + 1)) - set(cen.yr))
assert MISSING == [1998], MISSING   # decision year 1998 needs a 1997 evidence snapshot, which was never built
# First decision year plotted: 2000 (earlier years have small comparison sets; 1996, 1997 and 1999 are in the SI tables).
START = 2000
cen = cen[cen.yr >= START].copy(); Ty = cen[["yr"]].reset_index(drop=True)
XL = START - 1.0                     # left x-limit of the SI panels and the count strip
# absolute median percentiles (held-out pool): must equal point_estimates_test.tsv; pooled values are quoted in the SI
AB = pd.read_csv(f"{FR}/abs_percentiles.tsv", sep="\t")
AB = AB[(AB.pool == "test") & (AB.analysis == "pct_p1")]; AB["scope"] = AB.scope.astype(str)
PE = pd.read_csv(os.path.join(ROOT, "data", "derived", "si_event_time", "point_estimates_test.tsv"), sep="\t"); PE["scope"] = PE.scope.astype(str)
PE = PE[PE.analysis == "pct_p1"].set_index(["scope", "method"]).E1
def absser(m):
    a = AB[AB.method == m].set_index("scope")
    s = a.loc[[str(y) for y in Ty.yr]].reset_index(); s["yr"] = s.scope.astype(int)
    assert np.allclose(s.E1, [PE[(str(y), m)] for y in Ty.yr]), m
    return s, a.loc["pooled"]
ABS_D, POOL_D = absser("disatlas_ens5"); ABS_A, POOL_A = absser("adamic_adar_w")
assert round(POOL_D.E1, 2) == 24.00 and round(POOL_A.E1, 2) == 31.03, (POOL_D.E1, POOL_A.E1)

# ---- panel b: per-drug probability (AUC) that DisAtlas scores the next trialed disease above a
# popularity- and connectivity-matched alternative not yet trialed with the drug; held-out pool. Exact per event from the
# et02 percentile and matched-set size (trace_auc.py); yearly mean over drugs with 2,000-replicate drug bootstrap CIs.
AU = pd.read_csv(f"{FR}/auc_per_drug.tsv", sep="\t"); AU["scope"] = AU.scope.astype(str)
AU = AU[(AU.pool == "test") & (AU.method == "disatlas_ens5")]
def yser(an):
    Z = AU[AU.analysis == an].set_index("scope"); Z = Z.loc[[str(y) for y in Ty.yr]].reset_index(); Z["yr"] = Z.scope.astype(int); return Z
SETS = {"pct_p1": yser("pct_p1"), "pct_s3": yser("pct_s3"), "pct_s4": yser("pct_s4")}
AUY = SETS["pct_p1"]
CL = "#8FA9CF"   # lighter tint of the DisAtlas blue for the unmatched reference line
# Popularity-matched (solid, main colour) and no matching (lighter, dashed); held-out diseases.
PANEL_SER = [(SETS["pct_s4"], CL, "No matching", (0, (3, 1.6))), (SETS["pct_s3"], CD, "Popularity-matched", "-")]
for an, Z in SETS.items():
    print("panel b AUC %s: yearly %.3f-%.3f; years with CI above 0.5: %d of %d" % (an, Z.auc.min(), Z.auc.max(), int((Z.lo > 0.5).sum()), len(Z)))
DY = pd.read_csv(f"{FR}/auc_check_yearly_box.tsv", sep="\t")
for an in ("pct_p1", "pct_s3"):   # two scripts agree on yearly means and drug counts
    chk = DY[(DY.pool == "test") & (DY.analysis == an)].set_index("y").reindex(SETS[an].yr)
    assert np.allclose(chk["mean"].values, SETS[an].auc.values) and (chk["n_drugs"].values == SETS[an].n_drugs.values).all(), an
# ---- SI hits10 figure: pre-specified secondary endpoint, the share of
# drugs whose best (minimum) event percentile in the decision year is <= 10 (e2_per_drug.tsv); yearly mean
# over drugs with 2,000-replicate drug bootstrap CIs; held-out pool. Chance per drug = 1 - prod_i(1 - floor(N_i/10)/N_i) over
# its new trials of the year (exact for the discrete matched sets), averaged over the same drugs; drawn per year and set.
HT = pd.read_csv(f"{FR}/e2_per_drug_from2000.tsv", sep="\t"); HT["scope"] = HT.scope.astype(str)
HT = HT[(HT.pool == "test") & (HT.method == "disatlas_ens5") & (HT.k == 10)]
def hser(an):
    Z = HT[HT.analysis == an]
    H = Z[Z.kind == "e2"].set_index("scope").loc[[str(y) for y in Ty.yr]].reset_index()
    C = Z[Z.kind == "chc"].set_index("scope").loc[[str(y) for y in Ty.yr]].reset_index()
    H["yr"] = H.scope.astype(int); H["chance"] = C.value.values
    for c in ("value", "lo", "hi", "chance"): H[c] = 100 * H[c]
    assert (H.n.values == C.n.values).all(); return H
HSETS = {an: hser(an) for an in ("pct_p1", "pct_s3", "pct_s4")}
for an, Z in HSETS.items():
    print("panel b E2 %s: yearly %.1f-%.1f%%; chance %.1f-%.1f%%; years with CI above chance: %d of %d"
          % (an, Z.value.min(), Z.value.max(), Z.chance.min(), Z.chance.max(), int((Z.lo > Z.chance).sum()), len(Z)))
PANEL_HSER = [(HSETS["pct_s4"], CL, "No matching"), (HSETS["pct_s3"], CD, "Popularity-matched")]
WSER = None
if WEIGHTED:   # popularity-weighted single line
    WF = f"{FR}/weighted_yearly.tsv"
    assert WV in ("wex", "w51") and WU in ("event", "drug_year"), (WV, WU)
    COL = "hits10_event" if WU == "event" else "e2_drug_year"
    W0 = pd.read_csv(WF, sep="\t"); W0 = W0[(W0.pool == "test") & (W0.weights == WV) & (W0.trimmed == False)]
    def wser(m):
        W = W0[W0.method == m].rename(columns={"year": "yr"}).set_index("yr").loc[list(Ty.yr)].reset_index()
        return pd.DataFrame({"yr": W.yr, "value": 100 * W[COL], "lo": 100 * W[COL + "_lo"], "hi": 100 * W[COL + "_hi"],
                             "chance": 100 * W[COL + "_chance"], "n": W.n_drugs, "k": W.n_events})
    WSER, WPOP = wser("disatlas_ens5"), wser("popularity")
    assert (WSER.n.values == HSETS["pct_s3"].n.values).all(), "weighted file drug counts differ from the E2 trace"
    assert np.allclose(WSER.chance.values, WPOP.chance.values), "chance differs between methods"
    for nm, Z in (("DisAtlas", WSER), ("popularity", WPOP)):
        print("panel b weighted %s %s (%s, untrimmed): %s yearly %.1f-%.1f%%; chance %.1f-%.1f%%; CI above chance %d of %d; CI covers chance %d"
              % (WU, WV, "weighted_yearly.tsv", nm, Z.value.min(), Z.value.max(), Z.chance.min(), Z.chance.max(), int((Z.lo > Z.chance).sum()), len(Z),
                 int(((Z.lo <= Z.chance) & (Z.hi >= Z.chance)).sum())))
AUP = AU[(AU.analysis == "pct_p1") & (AU.scope == "pooled_drug_year")].iloc[0]; NABOVE = int((AUY.lo > 0.5).sum())
assert (AUY.n_drugs.values == HSETS["pct_s3"].n.values).all(), "P1 drug counts differ from the count strip (E2 trace)"
# PANEL   # p1 | hits (SI sensitivity: top-10% share, matched vs unmatched)
assert PANEL in ("p1", "hits", "or"), PANEL   # or = conditional-logit odds-ratio design
print("panel b AUC P1: pooled per drug-year %.3f [%.3f, %.3f]; years with CI above 0.5: %d of %d; yearly %.3f-%.3f"
      % (AUP.auc, AUP.lo, AUP.hi, NABOVE, len(AUY), AUY.auc.min(), AUY.auc.max()))

# ---- panel a data
TL = pd.read_csv(f"{NMO}/nmosd_timeline.tsv", sep="\t"); TL = TL[(TL.kind == "milestone") & (TL.drug == "eculizumab")].set_index("label")
FT = frac_year(TL.loc["First NMOSD trial", "date_or_slice"]); PA = frac_year(TL.loc["PNH FDA approval (BLA125166 ORIG-1)", "date_or_slice"])
NA = frac_year(TL.loc["NMOSD FDA approval (SUPPL-431)", "date_or_slice"])
# panel-a label years = calendar year of the full dates that place the markers; caption: March 2007, 2009, June 2019
for lab, x in [("PNH FDA approval (BLA125166 ORIG-1)", PA), ("First NMOSD trial", FT), ("NMOSD FDA approval (SUPPL-431)", NA)]:
    d = TL.loc[lab, "date_or_slice"]; assert len(d) == 10 and int(np.floor(x)) == pd.Timestamp(d).year, (lab, d, x)
assert [int(np.floor(x)) for x in (PA, FT, NA)] == [2007, 2009, 2019]
# literature milestones: plotted constants must agree with the sourced rows (tolerance 0.05 year)
LIT = pd.read_csv(f"{NMO}/nmosd_timeline.tsv", sep="\t", keep_default_na=False, dtype=str)
LIT = LIT[LIT.kind == "milestone"].set_index("label")
for lab, x in [("AQP4 antibody (NMO-IgG) splits NMO from MS", 2004.95), ("NMO criteria (2006)", 2006.4),
               ("NMOSD criteria (2015)", 2015.5), ("PNH pilot study (treatment start, approximate)", 2002.0)]:
    assert LIT.loc[lab, "ref_doi"] and abs(frac_year(LIT.loc[lab, "date_or_slice"]) - x) < 0.05, (lab, x)
assert pd.Timestamp(LIT.loc["PNH Phase 3 start (TRIUMPH)", "date_or_slice"]).year == 2004 and LIT.loc["PNH Phase 3 start (TRIUMPH)", "ref_doi"]
print("literature milestones sourced:", LIT.loc[LIT.ref_doi != "", "ref_citekey"].tolist())
Y = pd.read_csv(YEARLY, sep="\t"); Y = Y[(Y.pool == "test") & (Y.method == "disatlas_ens5")]
YA = pd.read_csv(YALL, sep="\t"); YA = YA[YA.pool == "test"]
_chk = YA[YA.method == "disatlas_ens5"].merge(Y, on=["y", "target"], suffixes=("", "_f"))
assert len(_chk) == len(Y) and np.allclose(_chk.pct_p1, _chk.pct_p1_f), "yearly tables disagree"
def rank_of(tgt, y, method="disatlas_ens5"):
    r = YA[(YA.method == method) & (YA.target == tgt) & (YA.y == y)].iloc[0]; n = int(r.n_p1); k = r.pct_p1 * n / 100
    assert abs(k - round(k)) < 1e-6; return int(round(k)), n
NY, MSY = 2007, 2002          # decision years drawn in the approved draft
RK, N = rank_of("NMO", NY); MSR, MSN = rank_of("MS", MSY)
AAR, AAN = rank_of("NMO", NY, "adamic_adar_w"); assert (AAR, AAN) == (3, 25), (AAR, AAN)   # illustrative case only
assert (RK, N) == (2, 25) and MSR <= 5  # draft wording "top 5" kept; exact rank printed below, (RK, N, MSR, MSN)
CASE = dict(N=N, RK=RK, aa_rank=AAR, aa_N=AAN, ms_rank=MSR, ms_N=MSN, ms_y=MSY, nmo_y=NY, first_trial=FT, pnh_approval=PA, nmosd_approval=NA)
# pool of N diseases includes NMO itself: the others are N - 1, of which N - RK rank below NMO
CASE["n_other"] = N - 1; CASE["n_lower"] = N - RK; CASE["pct_lower_of_others"] = 100 * (N - RK) / (N - 1)
assert (CASE["n_other"], CASE["n_lower"], round(CASE["pct_lower_of_others"])) == (24, 23, 96), CASE
print("panel a baseline: weighted Adamic-Adar NMO rank", AAR, "of", AAN, "decision year", NY, "(DisAtlas", RK, "of", N, ")")
assert MSN == 25 and MSR == 3, (MSR, MSN)   # MS: 3rd of 25, so "top 5" holds
print("panel a", CASE, "| NMO pct_p1 in table", float(Y[(Y.target == "NMO") & (Y.y == NY)].pct_p1.iloc[0]),
      "| step 2: %d like NMO + NMO = %d dots | step 3: %d of %d others ranked lower = %.1f%% | MS rank %d of %d (%d)"
      % (N - 1, N, N - RK, N - 1, CASE["pct_lower_of_others"], MSR, MSN, MSY))



def draw_panel_a_or(ax, CASE):
    """Panel a for the odds-ratio design: same timeline; ranks are among ALL not-yet-trialed held-out diseases (no matching)."""
    N, RK, MSR, MSN, MSY, NY = CASE["N"], CASE["RK"], CASE["ms_rank"], CASE["ms_N"], CASE["ms_y"], CASE["nmo_y"]
    FT, PA, NA = CASE["first_trial"], CASE["pnh_approval"], CASE["nmosd_approval"]
    ax.set_xlim(0,100); ax.set_ylim(0,57); ax.axis('off')
    X=lambda yr: 3+(yr-2000)/(2021-2000)*88
    y0=39; PU='#7B4F9E'; MSC='#5B3A80'
    for a_,b_,lab,fc in [(2000,2006.4,'treated as MS','#E4D9EE'),(2006.4,2015.5,'NMO','#CDB8E0'),(2015.5,2021,'NMOSD','#B497D1')]:
        ax.add_patch(Rectangle((X(a_),y0+7),X(b_)-X(a_),2.6,fc=fc,ec='white',lw=.8))
        ax.text((X(a_)+X(b_))/2,y0+8.3,lab,ha='center',va='center',fontsize=5.5,color='#3D2457',fontweight='bold')
    for yr,lab,ha,dx in [(2004.95,'AQP4 antibody\nsplits NMO from MS','right',.6),(2006.4,'NMO criteria\n(2006)','left',-.3),(2015.5,'NMOSD criteria\n(2015)','left',-.3)]:
        ax.plot([X(yr)]*2,[y0+9.6,y0+11],color=PU,lw=.6)
        ax.text(X(yr)+dx,y0+11.2,lab,ha=ha,va='bottom',fontsize=4.8,color=PU)
    # MS proxy bonus
    # drug timeline
    ax.annotate('',xy=(X(2021),y0),xytext=(X(2000),y0),arrowprops=dict(arrowstyle='->',lw=.9,color='k'))
    for yr in range(2000,2021,5): ax.text(X(yr),y0-1.3,str(yr),ha='center',va='top',fontsize=5.5,color='#555555'); ax.plot([X(yr)]*2,[y0-.5,y0+.5],color='k',lw=.6)
    # each label carries the year of the date that places its marker (exact dates, so no '~'; the pilot is approximate)
    yl = lambda x: int(np.floor(x))
    TL = []; ax._timeline_labels = TL   # checked for overlap in build_fig7_or
    for yr,lab,col in [(2002.0,'eculizumab\nPNH pilot (~2002)',ANC),(PA,f'eculizumab approved\nfor PNH ({yl(PA)})',ANC),
                       (FT,f'first NMOSD\ntrial ({yl(FT)})',CD),(NA,f'approved for\nNMOSD ({yl(NA)})',CD)]:
        ax.plot(X(yr),y0,'o',ms=3.5,color=col,zorder=4); ax.plot([X(yr)]*2,[y0,y0+1.3],color=col,lw=.6)
        TL.append(ax.text(X(yr)+{FT:-1.5,PA:-2.6}.get(yr,0),y0+1.4,lab,ha={FT:'left'}.get(yr,'center'),va='bottom',fontsize=5.5,color=col,linespacing=1.0))
    # DisAtlas markers below the axis: two brackets converging on the first trial
    xe=X(FT); FADE='#6A3D9A'; ax._bracket_labels = []   # saturated purple, classification-bar family
    YN = y0-4.6   # NMO bracket row
    for xs,yb_,col,lab,w,fw,ms in [
                                    (X(NY),YN,CD,f'≈{round(FT-NY)} years before first trial: NMO ranks {RK:,} of {N:,} by DisAtlas for eculizumab (top {CASE["top"]:.1f}%)',1.1,'bold',5.5)]:
        ax.plot(xs,yb_,marker='^',ms=ms,color=col if col!=CD else 'k',zorder=5)
        ax.plot([xs+.9,xe,xe],[yb_,yb_,y0-.3],color=col,lw=w)
        ax._bracket_labels.append(ax.text(xe+1,yb_+(0.9 if '\n' in lab else 0),lab,ha='left',va='center' if '\n' not in lab else 'top',fontsize=5.4,color=col,fontweight=fw))
    xp=X(NY); ym=YN
    yb=26.2
    ax.annotate('',xy=(xp,yb-.8),xytext=(xp,ym-1.6),arrowprops=dict(arrowstyle='->',lw=.8,color='#555555'))
    ax.text(xp+1,(yb+ym)/2-.8,'how this rank is computed',fontsize=5.2,va='center',ha='left',color='#555555',style='italic')
    # three-step strip
    def box(x,w,title,body,fc,ec):
        ax.add_patch(FancyBboxPatch((x,2),w,22,boxstyle='round,pad=0.3,rounding_size=1.5',fc=fc,ec=ec,lw=.8))
        ax.text(x+1.5,22.2,title,fontsize=6.5,fontweight='bold',va='top')
        return ax.text(x+1.5,4,body,fontsize=5.5,color='#444444',va='bottom')
    box(2,20,'1  Known use','PNH: pilot ~2002, phase 3\n2004, approved 2007','#FDF3EB',ANC)
    # antibody icon (IgG Y-shape)
    cx, cy = 7, 13
    for (x0_,y0_,x1_,y1_) in [(cx,cy-3.4,cx,cy),(cx,cy,cx-2.6,cy+3.6),(cx,cy,cx+2.6,cy+3.6)]:
        ax.plot([x0_,x1_],[y0_,y1_],color=ANC,lw=3.2,solid_capstyle='round')
        ax.plot([x0_,x1_],[y0_,y1_],color='#FDE9D9',lw=1.2,solid_capstyle='round')
    ax.text(11,16.6,'eculizumab',fontsize=6.3,fontweight='bold',va='center',color='k')
    ax.text(11,12.3,'antibody that\nblocks C5\n(complement)',fontsize=5,va='center',color='#444444')
    box(27,39,'2  Allow for popularity','NMO is compared with every disease not yet trialed\nwith eculizumab, and popularity is adjusted for','#F4F7FB','#888888')
    SC = CASE["scatter"]; rest = SC[~SC.is_nmo]; me = SC[SC.is_nmo].iloc[0]
    ins = ax.inset_axes([37.5, 11.6, 26.0, 6.6], transform=ax.transData)
    ins.scatter(rest.log1p_pop, rest.disatlas, s=1.2, c='#BBBBBB', lw=0, rasterized=True)
    ins.scatter([me.log1p_pop], [me.disatlas], s=10, c=CD, lw=0, zorder=4)
    ins.annotate('NMO', (me.log1p_pop, me.disatlas), xytext=(4, -1), textcoords='offset points', fontsize=5.2, color=CD, fontweight='bold', va='center')
    ins.set_xlabel('trials before 2007 (log scale)', fontsize=4.8, labelpad=1); ins.set_ylabel('DisAtlas\nscore', fontsize=4.8, labelpad=1)
    tv = [0, 1, 5, 20, 100]; ins.set_xticks(np.log1p(tv)); ins.set_xticklabels([str(v) for v in tv])
    ins.tick_params(labelsize=4.5, length=1.5, pad=1); ins.spines[['top', 'right']].set_visible(False)
    for sp in ins.spines.values(): sp.set_linewidth(.5)
    ins.set_facecolor('none'); ax._inset = ins
    ins.set_title(f'{N:,} candidate diseases, {int(NY)-1}', fontsize=5.0, color='#555555', pad=1.5)
    b3body = box(70,28,'3  Rank by DisAtlas score','score = closeness to PNH in the map\nno direct link to PNH in the map (2 hops apart)','#F4F7FB',CD)
    ax.text(71.5,17.6,f'NMO: rank {RK:,} of {N:,}',fontsize=6.5,color=CD,fontweight='bold',va='center')
    ax.text(71.5,15.2,f'= top {CASE["top"]:.1f}% of candidates',fontsize=5.5,va='center')
    # percentile strip, 0 = top; DisAtlas marker above, popularity alone below
    xs0, xs1, ys = 72.0, 95.5, 12.0
    ax.add_patch(Rectangle((xs0, ys - .45), xs1 - xs0, .9, fc='#E6E6E6', ec='none'))
    px = lambda pct: xs0 + (xs1 - xs0) * pct / 100
    ax.plot([px(50)] * 2, [ys - .45, ys + .45], color='#888888', lw=.8)
    B3 = [b3body]; ax._box3 = B3
    for xx, lab, ha in [(xs0, 'top', 'left'), (px(50), 'chance', 'center'), (xs1, 'bottom', 'right')]:
        B3.append(ax.text(xx, ys - 2.9, lab, fontsize=4.8, color='#777777', ha=ha, va='top'))
    ax.plot(px(CASE["top"]), ys + 1.15, marker='v', ms=4, color=CD, zorder=5)
    B3.append(ax.text(px(CASE["top"]) + 1.0, ys + 1.25, 'DisAtlas', fontsize=5.0, color=CD, va='center', fontweight='bold'))
    ax.plot(px(CASE["pop_top"]), ys - 1.15, marker='^', ms=3.4, color=G, zorder=5)
    B3.append(ax.text(px(CASE["pop_top"]) + 1.0, ys - 1.3, f'popularity alone: rank {CASE["pop_rank"]:,.0f} (top {CASE["pop_top"]:.0f}%)', fontsize=5.0, color='#666666', va='center'))
    for x0,x1 in [(22.6,26.3),(66.6,69.3)]: ax.annotate('',xy=(x1,13),xytext=(x0,13),arrowprops=dict(arrowstyle='->',lw=1,color='#555555'))



def draw_panel_b_p1(ax, Y):
    """Per-drug panel: DisAtlas per-drug probability that the drug's new trial disease ranks
    above a popularity- and connectivity-matched (P1) disease not yet trialed with it; yearly mean over drugs (dots joined
    by a line) with 95% drug-bootstrap CI bars (2,000 replicates); held-out pool; dashed 0.5 = random. DisAtlas only."""
    yy = np.arange(Y.yr.min(), Y.yr.max() + 1)
    ax.hlines(0.5, Y.yr.min() - .4, 2021.4, color='#6F6F6F', lw=.7, ls=(0, (4, 2)), zorder=1)
    ax.vlines(Y.yr, Y.lo, Y.hi, color=CD, lw=.8, alpha=.9, zorder=3)
    ax.plot(yy, Y.set_index('yr').auc.reindex(yy).values, color=CD, lw=.7, zorder=3); ax.plot(Y.yr, Y.auc, 'o', color=CD, ms=3.0, zorder=4)
    ax.text(2021.75, float(Y.auc.iloc[-1]), 'DisAtlas', fontsize=5.5, color=CD, va='center', ha='left')
    ax.text(2021.75, 0.5, 'random', fontsize=5.5, color='k', va='center', ha='left')
    lo, hi = 0.4, max(0.8, np.ceil(Y.hi.max() * 10 + 0.25) / 10)
    ax.set_ylim(lo, hi); ax.set_yticks(np.round(np.arange(lo, hi + 1e-9, .1), 1))
    ax.text(START, hi - 0.01, '\u2191 better', fontsize=5.5, color=GR, va='top', ha='left')
    ax.set_xlim(XL, 2024.5)
    ax.set_xlabel('Decision year (evidence up to the year before)')
    ax.set_ylabel('Probability\n(random = 0.5)', linespacing=1.0)
def draw_panel_b_weighted(ax, Y, P):
    """--weighted panel b: DisAtlas (dots, line, 95% CI bars; main colour), popularity alone under the same weights (grey;
    label 'popularity alone, same weights') and the chance line (thin dashed, labelled 'random'); direct labels at the right, no legend; 1998 never bridged."""
    yy = np.arange(Y.yr.min(), Y.yr.max() + 1)
    ax.plot(yy, Y.set_index('yr').chance.reindex(yy).values, color='#6F6F6F', lw=.7, ls=(0, (4, 2)), zorder=1)
    for Z, col, lw, ms, z in ((P, '#A8A8A8', .6, 2.4, 2), (Y, CD, .7, 3.0, 3)):
        ax.vlines(Z.yr, Z.lo, Z.hi, color=col, lw=.8, alpha=.9, zorder=z)
        ax.plot(yy, Z.set_index('yr').value.reindex(yy).values, color=col, lw=lw, zorder=z); ax.plot(Z.yr, Z.value, 'o', color=col, ms=ms, zorder=z + .5)
    top = np.ceil((Y.hi.max() + 5) / 10) * 10; ax.set_ylim(0, top); ax.set_yticks(range(0, int(top) + 1, 10 if top <= 60 else 20))
    L = sorted([(float(Y.value.iloc[-1]), 'DisAtlas', CD), (float(P.value.iloc[-1]), 'popularity alone,\nsame weights', '#7A7A7A'),
                (float(Y.chance.iloc[-1]), 'random', 'k')], key=lambda e: e[0])
    gap = 0.075 * top; pos = [L[0][0]]                    # centre-to-centre gap grows with the line count of both neighbours
    for p_, e in zip(L[:-1], L[1:]): pos.append(max(e[0], pos[-1] + gap * (1 + .6 * (p_[1].count('\n') + e[1].count('\n')))))
    for (v, t, col), yv in zip(L, pos): ax.text(2021.75, yv, t, fontsize=5.5, color=col, va='center', ha='left', linespacing=.95)
    ax.text(START, top * 0.95, '\u2191 better', fontsize=5.5, color=GR, va='top', ha='left')
    ax.set_xlim(XL, 2024.5)
    ax.set_xlabel('Decision year (evidence up to the year before)')
    ax.set_ylabel('New trials in the\ntop 10% (%)' if WU == 'event' else 'Drugs with a new trial\nin the top 10% (%)', linespacing=1.0)
def draw_panel_b_hits(ax, SER):
    """Hits10 panel: pre-specified secondary endpoint, the share of drugs whose best-ranked new trial of the year
    is in DisAtlas's top 10% of candidates, in percent; popularity-matched (solid, main colour) and no matching (lighter).
    Yearly mean over drugs as dots joined by a line, broken at 1998 (not scored: no model trained on data through 1997), 95%
    bootstrap CI as error bars, dots offset in x, direct labels at the right end. Chance is the random-ranking level for the
    same drugs, set and year (best of the drug's new trials; varies by year), thin dashed in the set's colour, labelled 'random'."""
    ends = []
    for k, (Y, col, lab) in enumerate(SER):
        dx = (-0.17, 0.17)[k]; x = Y.yr + dx
        yy = np.arange(Y.yr.min(), Y.yr.max() + 1)
        vv = Y.set_index('yr').value.reindex(yy).values; cc = Y.set_index('yr').chance.reindex(yy).values   # NaN at 1998
        z = 3 + 2 * k
        ax.plot(yy, cc, color=col if k else '#6F6F6F', lw=.7, ls=(0, (4, 2)), zorder=1)
        ax.vlines(x, Y.lo, Y.hi, color=col, lw=.8, alpha=.9, zorder=z)
        ax.plot(yy + dx, vv, color=col, lw=.7, zorder=z); ax.plot(x, Y.value, 'o', color=col, ms=3.0, zorder=z + 1)
        ends.append([float(Y.value.iloc[-1]), col, lab])
    ax.text(2021.75, (float(SER[0][0].chance.iloc[-1]) + float(SER[1][0].chance.iloc[-1])) / 2, 'random', fontsize=5.5, color='k', va='center', ha='left')   # one label for both chance lines
    for yv, col, lab in ends: ax.text(2021.75, yv, lab, fontsize=5.5, color=col, va='center', ha='left')
    ax.set_ylim(0, 90); ax.set_yticks(range(0, 91, 20))
    ax.text(1998 if START < 1998 else START, 85, '\u2191 better', fontsize=5.5, color=GR, va='top', ha='center' if START < 1998 else 'left')
    ax.set_xlim(XL, 2024.5)
    ax.set_xlabel('Decision year (evidence up to the year before)')
    ax.set_ylabel('Drugs with a new trial\nin the top 10% (%)', linespacing=1.0)
def draw_count_strip(ax, cen):
    ax.bar(cen.yr, cen.n_drugs, width=.7, color=CD, alpha=.35, lw=0)
    k=lambda n: f"{n/1000:.1f}k" if n>=1000 else f"{n}"
    for _,r in cen.iterrows():
        ax.text(r.yr, r.n_drugs+40, k(int(r.n_drugs)), ha='center', va='bottom', fontsize=4.3, color=CD)
        ax.text(r.yr, -120, k(int(r.n_events)), ha='center', va='top', fontsize=4.3, color=GR)
    ax.set_ylim(-520, max(1560, cen.n_drugs.max() + 160)); ax.axis('off')
    ax.text(2021.9, 700, 'unique drugs', ha='left', va='center', fontsize=5, color=CD)
    ax.text(2021.9, -330, 'drug–disease tests', ha='left', va='center', fontsize=5, color=GR)


# ---------------- odds-ratio design (post hoc conditional logit) ----------------
OR_LABEL = {"M2": "DisAtlas", "M2_gene_jaccard": "Shared genes", "M2_gobp_simgic": "GO processes", "M2_hpo_simgic": "HPO phenotypes",
            "M2_adamic_adar_w": "Weighted Adamic–Adar", "M2_strength": "Map connectivity"}
OR_TERM = {"M2": "z_dis", "M2_gene_jaccard": "z_gene_jaccard", "M2_gobp_simgic": "z_gobp_simgic", "M2_hpo_simgic": "z_hpo_simgic",
           "M2_adamic_adar_w": "z_adamic_adar_w", "M2_strength": "z_strength"}
def load_or(path):
    """Staged conditional-logit fits (Efron ties, full sample): one row per (set, model, term)."""
    A = pd.read_csv(path, sep="\t", dtype={"set": str})
    assert (A.groupby(["set", "model", "term"]).size() == 1).all(), "duplicate OR rows"
    return A, "full"
def org(A, model, term, st):
    r = A[(A.model == model) & (A.term == term) & (A.set == st)]; assert len(r) == 1, (model, term, st, len(r))
    r = r.iloc[0]; assert str(r.converged) in ("True", "TRUE", "1"), (model, st); return float(r.OR), float(r.lo95), float(r.hi95)
def draw_panel_b_or(ax, A, years):
    SP = len(years) <= 4; xs = list(range(len(years))) if SP else list(years); xp = (len(years) - 1 + 1.3) if SP else max(xs) + 3.0
    ser = [("M0", "z_dis", '#9DB3DA', 'o', 3.2, "DisAtlas, unadjusted", -0.35), ("M2", "z_dis", CD, 'o', 4.6, "DisAtlas, adjusted for popularity", 0.0),
           ("PERM_M2", "z_dis_perm", '#8A8A8A', 's', 3.2, "shuffled DisAtlas scores, adjusted", 0.35)]
    ax.axhline(1, color='k', lw=.7, ls=(0, (4, 3)), zorder=1)
    out = {}
    for m, tm, col, mk, ms, lab, dx in ser:
        v = [org(A, m, tm, str(y)) for y in years]; p = org(A, m, tm, "pooled"); out[m] = (v, p)
        dx_ = dx * (0.45 if SP else 1.0)
        ax.errorbar([x + dx_ for x in xs], [a for a, _, _ in v], yerr=[[a - l for a, l, _ in v], [h - a for a, _, h in v]], fmt=mk, ms=ms, color=col,
                    lw=1.0 if m == "M2" else .8, capsize=0, zorder=3 if m == "M2" else 2, label=lab)
        if len(xs) > 4 and m == "M2": ax.plot(xs, [a for a, _, _ in v], color=col, lw=.6, zorder=2)
        ax.errorbar([xp + dx_], [p[0]], yerr=[[p[0] - p[1]], [p[2] - p[0]]], fmt=mk, ms=ms, color=col, lw=1.0, capsize=0, zorder=3)
    ax._labels = []
    v, p = out["M2"]
    if len(xs) <= 4:
        for x, (a, _, _) in zip(xs, v): ax._labels.append(ax.text(x + 0.12, a, f'{a:.2f}', fontsize=6, color=CD, fontweight='bold', va='center'))
    ax._labels.append(ax.text(xp + (0.12 if SP else 0.55), p[0], f'{p[0]:.2f}', fontsize=6, color=CD, fontweight='bold', va='center'))
    ax.axvline((max(xs) + xp) / 2, color='#CCCCCC', lw=.6)
    ax._noeff = ax.text(1.012, 1.0, 'no\neffect', fontsize=5.5, va='center', ha='left', color='#333333', linespacing=1.0, clip_on=False,
                        transform=mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData))
    allv = [x for m in out for (a, l, h) in out[m][0] + [out[m][1]] for x in (l, h)]
    ax.set_yscale('log'); ax.set_ylim(min(0.8, min(allv) * 0.95), max(7.0, max(allv) * 1.6))
    tk = [0.8, 1, 1.5, 2, 3, 4, 5]; ax.set_yticks(tk); ax.set_yticklabels([f'{t:g}' for t in tk]); ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    yt = list(xs) if SP else [y for y in xs if (y % 4 == 0 or y == min(xs))]; ytl = [str(years[k]) for k in yt] if SP else [str(y) for y in yt]
    pl = 'pooled' if len(xs) > 4 else f'{len(xs)} years\npooled'
    ax.set_xticks(yt + [xp]); ax.set_xticklabels(ytl + [pl]); ax.set_xlim(min(xs) - (0.55 if SP else 1.8), xp + (0.6 if SP else 2.0))
    ax.set_xlabel('Decision year (evidence up to the year before)'); ax.set_ylabel('Odds ratio per 1 SD of score\n(log scale)')
    ax.legend(loc='upper left', fontsize=5.5, frameon=False, handletextpad=0.3, borderaxespad=0.2)
    return out
def draw_panel_c_or(ax, A, npool):
    rows = [(m, *org(A, m, OR_TERM[m], "pooled")) for m in OR_LABEL]
    rows.sort(key=lambda r: -r[1]); assert rows[0][0] == "M2", "DisAtlas is not the largest pooled OR"
    ax.axvline(1, color='k', lw=.7, ls=(0, (4, 3)))
    ax._labels = []
    for i, (m, o, l, h) in enumerate(rows):
        yy = len(rows) - 1 - i; dis = m == "M2"; col = CD if dis else '#5A5A5A'
        ax._ci = getattr(ax, '_ci', []) + [(m, l, h, yy, 4.4 if dis else 3.2)]
        ms_ = 4.4 if dis else 3.2
        ax.plot([o], [yy], 'o', ms=ms_, color=col, mec='white', mew=.5, zorder=3)
        # 95% CI drawn ON TOP of the marker with caps taller than the marker, so even the narrowest interval is visible
        ax.errorbar([o], [yy], xerr=[[o - l], [h - o]], fmt='none', ecolor='k' if dis else '#222222', elinewidth=.9, capsize=ms_ * 0.95, capthick=.9, zorder=4)
        ax._labels.append(ax.text(h * 1.025, yy, f'{o:.2f}', fontsize=6 if dis else 5.6, color=col, fontweight='bold' if dis else 'normal', va='center'))
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([OR_LABEL[r[0]] for r in rows][::-1])
    for t in ax.get_yticklabels():
        if t.get_text() == "DisAtlas": t.set_color(CD); t.set_fontweight('bold')
    ax.set_xscale('log'); ax.set_xlim(min(0.92, min(r[2] for r in rows) * 0.97), max(2.25, max(r[3] for r in rows) * 1.12)); xt = [1, 1.25, 1.5, 2]; ax.set_xticks(xt); ax.set_xticklabels([f'{t:g}' for t in xt])
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter()); ax.set_ylim(-0.6, len(rows) - 0.4)
    ax.set_xlabel('Odds ratio per 1 SD, adjusted for popularity\n(' + npool + ')')
    return rows
def build_fig7_or():
    fig = plt.figure(figsize=(7.2, 6.4))
    axa = fig.add_axes([0.0, 0.47, 1.0, 0.53]); draw_panel_a_or(axa, CASE_OR)
    axb = fig.add_axes([0.085, 0.085, 0.37, 0.30]); axc = fig.add_axes([0.665, 0.085, 0.30, 0.30])
    yrs = sorted(int(v) for v in ORA.set.unique() if v != "pooled" and int(v) >= START)
    OUTB = draw_panel_b_or(axb, ORA, yrs); OUTC = draw_panel_c_or(axc, ORA, f'{len(yrs)} years pooled' if len(yrs) <= 4 else f'{min(yrs)}–{max(yrs)} pooled')
    for ax in (axb, axc):
        ax.spines[['top', 'right']].set_visible(False); ax.tick_params(labelsize=6); ax.xaxis.label.set_size(6.5); ax.yaxis.label.set_size(6.5)
    tb_ = axb.set_title("Odds that a disease is the drug's next trial", fontsize=6.8, loc='left'); tc_ = axc.set_title("Each method, same adjustment", fontsize=6.8, loc='left')
    fig.canvas.draw(); r_ = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    ca = [inv.transform(a.get_window_extent(r_)) for a in list(axa.texts) + list(axa.patches) + list(axa.lines)]
    x0a = min(b[0, 0] for b in ca); y0a = min(b[0, 1] for b in ca); y1a = max(b[1, 1] for b in ca)
    fig.text(x0a - 0.03, y1a, 'a', fontsize=9, fontweight='bold', va='top')
    bb = inv.transform(axb.get_tightbbox(r_)); fig.text(min(x0a - 0.03, bb[0, 0] - 0.005), bb[1, 1], 'b', fontsize=9, fontweight='bold', va='top')
    bc = inv.transform(axc.get_tightbbox(r_)); fig.text(bc[0, 0] - 0.01, bc[1, 1], 'c', fontsize=9, fontweight='bold', va='top')
    assert bb[1, 1] < y0a - 0.01 and bc[1, 1] < y0a - 0.01, ('panels b/c collide with panel a', bb[1, 1], bc[1, 1], y0a)
    assert bb[1, 0] < bc[0, 0] - 0.01, ('panel b runs into panel c', bb[1, 0], bc[0, 0])
    if ORSRC == "smoke":
        fig.canvas.draw(); r_ = fig.canvas.get_renderer(); ylo = min(inv.transform(ax.get_tightbbox(r_))[0, 1] for ax in (axb, axc))
        fig.text(0.5, ylo - 0.012, f"DRAFT from smoke test (decision years {', '.join(map(str, yrs))}; all not-yet-trialed held-out diseases; conditional logit, Efron). Full 2000–2021 run pending.",
                 fontsize=5.5, color='#B22222', style='italic', ha='center', va='top')
    fig.canvas.draw(); r_ = fig.canvas.get_renderer()
    TL = [t.get_window_extent(r_) for t in axa._timeline_labels]
    for i in range(len(TL)):
        for j in range(i + 1, len(TL)): assert not TL[i].overlaps(TL[j]), ('timeline labels overlap', i, j)
    BL = [t.get_window_extent(r_) for t in axa._bracket_labels]; FW = fig.bbox
    assert len(BL) == 1, 'expected only the NMO bracket label (MS line removed)'
    # panel c: every 95% CI narrower than its marker -> state it in the panel (caption carries the same text)
    dpi = fig.dpi; wci = [(m, abs(axc.transData.transform((h, yy))[0] - axc.transData.transform((l, yy))[0]), ms * dpi / 72) for m, l, h, yy, ms in axc._ci]
    print('panel c CI width vs marker (px):', [(m, round(w, 1), round(d, 1)) for m, w, d in wci])
    # CIs are drawn visibly (capped, above markers); caps exceed the marker height so every interval is legible
    fig.canvas.draw(); r_ = fig.canvas.get_renderer()
    # 'no effect' label: clear of every marker/errorbar in panel b and of the legend
    ne = axb._noeff.get_window_extent(r_)
    for coll in axb.lines + axb.collections:
        try: pts = coll.get_xydata() if hasattr(coll, 'get_xydata') else None
        except Exception: pts = None
        if pts is None or not len(pts) or coll.get_linestyle() in ('--',) : continue
        P = axb.transData.transform(pts)
        if coll.get_marker() not in (None, 'None', '', ' '):
            assert not any(ne.x0 - 4 < x < ne.x1 + 4 and ne.y0 - 4 < y < ne.y1 + 4 for x, y in P), ("'no effect' touches a marker", coll.get_label())
    for c_ in axb.collections:   # errorbar line collections
        for seg in c_.get_segments():
            P = axb.transData.transform(seg)
            if len(P) == 2:
                (xa, ya), (xb, yb_) = P
                if ne.x0 - 2 < xa < ne.x1 + 2 and min(ya, yb_) < ne.y1 and max(ya, yb_) > ne.y0: raise AssertionError("'no effect' touches an error bar")
    assert not ne.overlaps(axb.get_legend().get_window_extent(r_)), "'no effect' overlaps legend"
    BL = [t.get_window_extent(r_) for t in axa._bracket_labels]
    for b_ in BL:
        assert FW.x0 < b_.x0 and b_.x1 < FW.x1, ('bracket label outside figure', b_.x0, b_.x1, FW.x1)
        for t in TL: assert not b_.overlaps(t), 'bracket label overlaps timeline label'
    print('bracket labels right edge (px):', [round(b_.x1) for b_ in BL], 'figure width', round(FW.x1))
    TA = [t for t in axa.texts]; BA = [t.get_window_extent(r_) for t in TA]
    insb = axa._inset.get_tightbbox(r_)
    for t, b_ in zip(TA, BA):
        if t.get_text() != 'NMO' and b_.overlaps(insb) and t.get_position()[1] < 20 and 27 < t.get_position()[0] < 66: raise AssertionError(('box-2 text overlaps scatter', t.get_text()))
    B3 = [t.get_window_extent(r_) for t in axa._box3]
    for i in range(len(B3)):
        for j in range(i + 1, len(B3)): assert not B3[i].overlaps(B3[j]), ('box-3 overlap', axa._box3[i].get_text(), axa._box3[j].get_text())
    for ax in (axb, axc):    # value labels vs each other, the legend and tick labels
        T = ax._labels + ([ax.get_legend()] if ax.get_legend() else []); B = [t.get_window_extent(r_) for t in T]
        for i in range(len(B)):
            for j in range(i + 1, len(B)): assert not B[i].overlaps(B[j]), ('label overlap', i, j)
    hb = [t.get_window_extent(r_) for t in axb.get_legend().get_texts()]
    for (v, p) in [OUTB["M0"]]: pass
    from matplotlib.transforms import Bbox
    bbs = [a.get_window_extent(r_) for a in list(axa.texts) + list(axa.patches) + list(axa.lines) + list(fig.texts)] + [ax.get_tightbbox(r_) for ax in (axb, axc, axa._inset)]
    U = Bbox.union(bbs).transformed(fig.dpi_scale_trans.inverted())
    fig._crop_bbox = Bbox.from_extents(U.x0 - .05, U.y0 - .05, U.x1 + .05, U.y1 + .05)
    print("panel b M2 yearly:", [(y, round(a, 3), round(l, 3), round(h, 3)) for y, (a, l, h) in zip(yrs, OUTB["M2"][0])], "pooled", tuple(round(x, 3) for x in OUTB["M2"][1]))
    print("panel b M0 pooled", tuple(round(x, 3) for x in OUTB["M0"][1]), "PERM_M2 pooled", tuple(round(x, 3) for x in OUTB["PERM_M2"][1]))
    print("panel c pooled:", [(OR_LABEL[m], round(o, 3), round(l, 3), round(h, 3)) for m, o, l, h in OUTC])
    return fig
if PANEL == "or":
    import json as _json
    ORA, ORSRC = load_or(f"{OD}/or_clogit.tsv")
    CJ = _json.load(open(f"{OD}/nmo_case_clogit.json")); SCT = pd.read_csv(f"{OD}/nmo_case_scatter.tsv", sep="\t")
    assert CJ["year"] == CASE["nmo_y"] and CJ["nmo_in_test_pool"] and CJ["disatlas_ties"] == 0 and len(SCT) == CJ["N"] == CJ["yearly_table_n_cand"]
    assert abs(CJ["yearly_table_pct_s4"] - CJ["disatlas_pct_rank"]) < 1e-9 and SCT.is_nmo.sum() == 1
    msr = YA[(YA.method == "disatlas_ens5") & (YA.target == "MS") & (YA.y == CASE["ms_y"])].iloc[0]
    MSK = msr.pct_s4 * msr.n_cand / 100; assert abs(MSK - round(MSK)) < 1e-6, MSK
    CASE_OR = dict(CASE, N=int(CJ["N"]), RK=int(CJ["disatlas_rank"]), top=CJ["disatlas_pct_rank"], pop_rank=CJ["popularity_rank"],
                   pop_top=CJ["popularity_pct_rank"], ms_rank=int(round(MSK)), ms_N=int(msr.n_cand), ms_top=float(msr.pct_s4), scatter=SCT)
    print("panel a (odds-ratio design): NMO", CASE_OR["RK"], "of", CASE_OR["N"], "top %.2f%%" % CASE_OR["top"], "| popularity alone rank",
          CASE_OR["pop_rank"], "top %.2f%%" % CASE_OR["pop_top"], "| MS", CASE_OR["ms_rank"], "of", CASE_OR["ms_N"], "| OR source", ORSRC)

def build_panel_b_only(cen):
    fig = plt.figure(figsize=(7.2, 2.75))
    axs = fig.add_axes([0.11, 0.70, 0.74, 0.17]); axb = fig.add_axes([0.11, 0.14, 0.74, 0.52], sharex=axs)
    draw_panel_b_weighted(axb, WSER, WPOP) if WEIGHTED else (draw_panel_b_p1(axb, AUY) if PANEL == 'p1' else draw_panel_b_hits(axb, PANEL_HSER))
    axb.spines[['top', 'right']].set_visible(False); axb.tick_params(labelsize=6); axb.xaxis.label.set_size(6.5); axb.yaxis.label.set_size(6.5)
    axb.set_xticks(range(START + (START % 2), 2022, 2)); axb.set_xlim(XL + 0.2, 2024.5)
    draw_count_strip(axs, cen)
    fig.canvas.draw(); r_ = fig.canvas.get_renderer()
    from matplotlib.transforms import Bbox
    bbs = [ax.get_tightbbox(r_) for ax in (axs, axb)]
    U = Bbox.union(bbs).transformed(fig.dpi_scale_trans.inverted())
    fig._crop_bbox = Bbox.from_extents(U.x0 - .05, U.y0 - .05, U.x1 + .05, U.y1 + .05)
    T = [t for ax in (axs, axb) for t in ax.texts]; bb = [t.get_window_extent(r_) for t in T]
    for i in range(len(bb)):                                   # render-then-verify: no text overlaps
        for j in range(i + 1, len(bb)): assert not bb[i].overlaps(bb[j]), ('text overlap', T[i].get_text(), T[j].get_text())
    return fig
fig = build_fig7_or() if PANEL == "or" else build_panel_b_only(cen)
os.makedirs(OUT, exist_ok=True)
STEM = {None: "fig7_prospective", "hits10": "SI_event_time_hits10", "weighted": "SI_event_time_weighted"}[SI]
fig.savefig(f"{OUT}/{STEM}.png", dpi=300, bbox_inches=fig._crop_bbox, facecolor="white")
fig.savefig(f"{OUT}/{STEM}.pdf", bbox_inches=fig._crop_bbox, facecolor="white")
if SI is None:   # value table read by tests/test_main_parity.py (keys fig7_prospective:<key>); one row per main-text number
    import json as _json
    _o = pd.read_csv(os.path.join(OD, "or_clogit.tsv"), sep="\t", dtype={"set": str}); _o = _o[_o.model != "H2H"].set_index(["set", "model"], verify_integrity=True)
    _V = []
    def _or(key, s, m, ci=True):
        r = _o.loc[(s, m)]; _V.append((key, r.OR))
        if ci: _V.extend([(key.replace("|", "_lo|") if "|" in key else key + "_lo", r.lo95), (key.replace("|", "_hi|") if "|" in key else key + "_hi", r.hi95)])
        return r
    _p = _or("or_disatlas", "pooled", "M2"); _V += [("n_trials_or", _p.n_cases), ("n_drugs_or", _p.n_drugs)]
    assert all(_o.loc[("pooled", m)].n_cases == _p.n_cases for m in ("PERM_M2", "M2_strength", "M2_gene_jaccard", "M2_adamic_adar_w"))
    _or("or_shuffled", "pooled", "PERM_M2"); _or("or_connectivity", "pooled", "M2_strength")
    _or("or_gene", "pooled", "M2_gene_jaccard"); _or("or_waa", "pooled", "M2_adamic_adar_w")
    _or("or_gobp", "pooled", "M2_gobp_simgic", ci=False); _or("or_hpo", "pooled", "M2_hpo_simgic", ci=False)
    _y = {int(s): _o.loc[(s, "M2")].OR for s, m in _o.index if m == "M2" and s.isdigit() and int(s) >= START}
    assert min(_y, key=_y.get) == 2020                      # 'It is lowest in 2020'
    _or("or_yearly|2020", "2020", "M2"); _or("or_yearly|2021", "2021", "M2", ci=False)
    _c = _json.load(open(os.path.join(OD, "nmo_case_clogit.json")))
    _h = _json.load(open(os.path.join(OD, "nmo_hops.json")))["stripped_training_graph"]; assert _h["direct_edge"] is False
    _cs = pd.read_csv(os.path.join(OD, "nmo_candidate_sets.tsv"), sep="\t")
    _al = _cs[(_cs.pool != "held-out pool") & (_cs.set == "any distance (current)")].iloc[0]
    _V += [("nmo_rank_disatlas", _c["disatlas_rank"]), ("nmo_pool", _c["N"]), ("nmo_top_pct_disatlas", _c["disatlas_pct_rank"]),
           ("nmo_rank_popularity", _c["popularity_rank"]), ("nmo_top_pct_popularity", _c["popularity_pct_rank"]), ("nmo_hops_pnh", _h["hops"]),
           ("nmo_rank_all", _al.disatlas_rank), ("nmo_pool_all", _al.N), ("nmo_top_pct_all", 100 * _al.disatlas_rank / _al.N)]
    assert len({k for k, _ in _V}) == len(_V) == 32
    pd.DataFrame(_V, columns=["key", "value"]).to_csv(f"{OUT}/fig7_prospective_values.tsv", sep="\t", index=False)
print("wrote", OUT)
