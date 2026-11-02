"""Main Figure 7: simulated prospective (event-time) validation.
a, post hoc case: eculizumab, PNH -> NMO(SD); b, all-drug yearly advantage of DisAtlas over weighted Adamic-Adar
(positive = DisAtlas better) with per-year drug / test counts.
Inputs (data/derived/fig7_prospective/):
  per_year_vs_wAA.tsv, deltas.tsv, merge_and_matching_census.tsv   event-time results (test pool)
  nmosd_timeline.tsv                                               dated milestones with sources
  fig7_nmo_ms_yearly_ens5_test.tsv                                 five-seed-mean NMO / MS yearly ranks (test pool)
  nmo_ms_yearly_all_methods.tsv                                    all methods, all pools (panel a baseline rank)
  milestone_reference_verification.tsv                             reference check for the literature dates
Literature dates (PNH pilot ~2002, phase 3 2004, AQP4 antibody 2004, NMO criteria 2006, NMOSD criteria 2015) are drawn
as constants and asserted against the sourced milestone rows of nmosd_timeline.tsv.
    python figures/main/fig7_prospective.py [OUTDIR]
Writes fig7_prospective.png/.pdf and fig7_prospective_values.tsv next to this script.
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
OUT = sys.argv[1] if len(sys.argv) > 1 else HERE
ET = NMO = os.path.join(ROOT, "data", "derived", "fig7_prospective")
YEARLY = os.path.join(ET, "fig7_nmo_ms_yearly_ens5_test.tsv")
YALL = os.path.join(ET, "nmo_ms_yearly_all_methods.tsv")   # all methods, all pools; baseline rank in panel a
CD='#1F4E99'; PUf='#A893C4'; GR='#555555'; G='#9E9E9E'; ANC='#B35806'
ordn = lambda k: {1: "1st", 2: "2nd", 3: "3rd"}.get(k, f"{k}th")
def frac_year(d):
    d = pd.Timestamp(d); return d.year + (d.dayofyear - 1) / 365.25

# ---- panel b data
PY = pd.read_csv(f"{ET}/per_year_vs_wAA.tsv", sep="\t")
PY = PY[(PY.pool == "test") & (PY.analysis == "pct_p1") & (PY.comparator == "adamic_adar_w")].copy()
PY = PY[PY.scope.astype(str).str.isdigit()]; PY["yr"] = PY.scope.astype(int); Ty = PY.sort_values("yr")
DL = pd.read_csv(f"{ET}/deltas.tsv", sep="\t")
pr = DL[(DL.pool == "test") & (DL.scope == "pooled") & (DL.analysis == "pct_p1") & (DL.comparator == "adamic_adar_w")].iloc[0]
pooled = dict(d=pr.dE1, lo=pr.dE1_lo, hi=pr.dE1_hi)
dly = DL[(DL.pool == "test") & (DL.analysis == "pct_p1") & (DL.comparator == "adamic_adar_w") & DL.scope.astype(str).str.isdigit()].copy()
dly["yr"] = dly.scope.astype(int); dly = dly.set_index("yr").loc[Ty.yr]
assert np.allclose(dly[["dE1", "dE1_lo", "dE1_hi"]].to_numpy(), Ty[["dE1", "dE1_lo", "dE1_hi"]].to_numpy()), "per-year table != deltas"
CE = pd.read_csv(f"{ET}/merge_and_matching_census.tsv", sep="\t"); CE = CE[CE.pool == "test"]
NDRUG = int(CE[CE.scope == "pooled"].n_drugs.iloc[0])
cen = CE[CE.scope.astype(str).str.isdigit()].copy(); cen["yr"] = cen.scope.astype(int); cen = cen[cen.yr.isin(Ty.yr)].sort_values("yr")
assert cen.yr.tolist() == Ty.yr.tolist()
# manuscript values
assert NDRUG == 6496 and round(-pooled["d"], 1) == 7.0, (NDRUG, pooled)
MISSING = sorted(set(range(Ty.yr.min(), Ty.yr.max() + 1)) - set(Ty.yr))
print("panel b years", Ty.yr.min(), Ty.yr.max(), "n", len(Ty), "missing", MISSING, "pooled", {k: round(-v, 2) for k, v in pooled.items()},
      "drugs", NDRUG, "events", int(CE[CE.scope == "pooled"].n_events.iloc[0]))
NCI, NYRS = int((-Ty.dE1_hi > 0).sum()), len(Ty)
assert MISSING == [1998], MISSING   # decision year 1998 needs a 1997 evidence snapshot, which was never built
print("years with CI above 0:", int((-Ty.dE1_hi > 0).sum()), "below 0:", int((-Ty.dE1_lo < 0).sum()), "point>0:", int((-Ty.dE1 > 0).sum()))

# ---- panel a data
TL = pd.read_csv(f"{NMO}/nmosd_timeline.tsv", sep="\t"); TL = TL[(TL.kind == "milestone") & (TL.drug == "eculizumab")].set_index("label")
FT = frac_year(TL.loc["First NMOSD trial", "date_or_slice"]); PA = frac_year(TL.loc["PNH FDA approval (BLA125166 ORIG-1)", "date_or_slice"])
NA = frac_year(TL.loc["NMOSD FDA approval (SUPPL-431)", "date_or_slice"])
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
NY, MSY = 2007, 2002          # decision years drawn in panel a
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

def draw_panel_a(ax, CASE):
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
    for yr,lab,col in [(2002.0,'eculizumab\nPNH pilot (~2002)',ANC),(PA,'eculizumab\napproved for PNH',ANC),(FT,'first NMOSD\ntrial',CD),(NA,'approved\nfor NMOSD',CD)]:
        ax.plot(X(yr),y0,'o',ms=3.5,color=col,zorder=4); ax.plot([X(yr)]*2,[y0,y0+1.3],color=col,lw=.6)
        ax.text(X(yr)+{PA:1.5,FT:-1.5}.get(yr,0),y0+1.4,lab,ha={PA:'right',FT:'left'}.get(yr,'center'),va='bottom',fontsize=5.5,color=col)
    # DisAtlas markers below the axis: two brackets converging on the first trial
    xe=X(FT); FADE='#A893C4'
    for xs,yb_,col,lab,w,fw,ms in [(X(MSY),y0-4.0,FADE,f'≈{round(FT-MSY)} years before first trial: MS (proxy for NMO) ranks top 5 by DisAtlas for eculizumab',.8,'normal',4.5),
                                    (X(NY),y0-7.4,CD,f'≈{round(FT-NY)} years before first trial: NMO ranks {ordn(RK)} of {N} by DisAtlas for eculizumab',1.1,'bold',5.5)]:
        ax.plot(xs,yb_,marker='^',ms=ms,color=col if col!=CD else 'k',zorder=5)
        ax.plot([xs+.9,xe,xe],[yb_,yb_,y0-.3],color=col,lw=w)
        ax.text(xe+1,yb_+(0.9 if '\n' in lab else 0),lab,ha='left',va='center' if '\n' not in lab else 'top',fontsize=5.4 if fw=='bold' else 5.0,color=col,fontweight=fw)
    xp=X(NY); ym=y0-7.4
    yb=26.2
    ax.annotate('',xy=(xp,yb-.8),xytext=(xp,ym-1.6),arrowprops=dict(arrowstyle='->',lw=.8,color='#555555'))
    ax.text(xp+1,(yb+ym)/2-.8,'how this rank is computed',fontsize=5.2,va='center',ha='left',color='#555555',style='italic')
    # three-step strip
    def box(x,w,title,body,fc,ec):
        ax.add_patch(FancyBboxPatch((x,2),w,22,boxstyle='round,pad=0.3,rounding_size=1.5',fc=fc,ec=ec,lw=.8))
        ax.text(x+1.5,22.2,title,fontsize=6.5,fontweight='bold',va='top')
        ax.text(x+1.5,4,body,fontsize=5.5,color='#444444',va='bottom')
    box(2,20,'1  Known use','PNH: pilot ~2002, phase 3\n2004, approved 2007','#FDF3EB',ANC)
    # antibody icon (IgG Y-shape)
    cx, cy = 7, 13
    for (x0_,y0_,x1_,y1_) in [(cx,cy-3.4,cx,cy),(cx,cy,cx-2.6,cy+3.6),(cx,cy,cx+2.6,cy+3.6)]:
        ax.plot([x0_,x1_],[y0_,y1_],color=ANC,lw=3.2,solid_capstyle='round')
        ax.plot([x0_,x1_],[y0_,y1_],color='#FDE9D9',lw=1.2,solid_capstyle='round')
    ax.text(11,16.6,'eculizumab',fontsize=6.3,fontweight='bold',va='center',color='k')
    ax.text(11,12.3,'antibody that\nblocks C5\n(complement)',fontsize=5,va='center',color='#444444')
    box(27,39,'2  Control for popularity bias','busy diseases get new trials anyway,\nso they are removed before ranking','#F4F7FB','#888888')
    rng = np.random.default_rng(5)
    S = rng.choice([30, 16, 7, 3], size=(4, 6), p=[.25, .25, .25, .25])
    for i in range(4):
        for j in range(6):
            ax.scatter(30.5 + j * 2.6, 9.4 + i * 1.95, s=S[i, j], c='#BBBBBB', lw=0)
    ax.text(37, 16.6, f'all diseases, {int(NY)-1}\n(dot size = trials)', ha='center', fontsize=5.2, color='#555555')
    ax.annotate('', xy=(48.5, 12.3), xytext=(45.5, 12.3), arrowprops=dict(arrowstyle='->', lw=.8, color='#555555'))
    k = 0
    for i in range(4):
        for j in range(7):
            if k >= N: break
            nmo = (i == 1 and j == 3)
            ax.scatter(50.5 + j * 2.0, 9.4 + i * 1.95, s=9 if nmo else 4, c=CD if nmo else '#7F7F7F', lw=0, zorder=4); k += 1
    ax.text(56.5, 16.6, f'{CASE["n_other"]} like NMO: similar\ntrials and links', ha='center', fontsize=5.2, color='#555555')
    box(70,28,'3  Rank by DisAtlas score',f'score = closeness to PNH in the map\nweighted Adamic\u2013Adar: NMO {AAR} of {AAN}','#F4F7FB',CD)
    for r in range(1,6):
        yy=18.5-(r-1)*2.3; nmo=r==RK
        ax.text(74.6,yy,f'{r}',fontsize=5.5,va='center',ha='right',color='#666666')
        ax.add_patch(Rectangle((75.5,yy-.7),(9.5-1.2*r),1.4,color=CD if nmo else G,lw=0))
    ax.text(86.5,17.4,f'NMO: rank\n{RK} of {N}',fontsize=6.5,color=CD,fontweight='bold',va='center')
    ax.text(86.5,14.85,f'= {CASE["pct_lower_of_others"]:.0f}% of others\nranked lower\n(chance 50%)',fontsize=5.5,va='top',linespacing=1.05)
    for x0,x1 in [(22.6,26.3),(66.6,69.3)]: ax.annotate('',xy=(x1,13),xytext=(x0,13),arrowprops=dict(arrowstyle='->',lw=1,color='#555555'))


def draw_panel_c(ax, Ty, pooled):
    adv=-Ty.dE1; lo=-Ty.dE1_hi; hi=-Ty.dE1_lo
    ax.axhspan(-pooled['hi'],-pooled['lo'],color=CD,alpha=.12,lw=0)
    ax.axhline(-pooled['d'],color=CD,lw=.9)
    ax.text(2021.6,-pooled['lo']+0.4,f"all years\ncombined\n+{-pooled['d']:.1f}",fontsize=5.5,color=CD,va='bottom',ha='left')
    ax.vlines(Ty.yr,lo,hi,color=CD,lw=.8,alpha=.8); ax.plot(Ty.yr,adv,'o',color=CD,ms=3.2,zorder=3)
    ax.axhline(0,color='k',lw=.7)
    ax.text(2024.3,-0.8,'below 0:\nAdamic–Adar\nbetter',fontsize=5.5,color=GR,va='top',ha='right')
    ax.set_xlim(1995,2024.5); ax.set_xticks(range(1996,2022,4))
    ax.set_xlabel('Decision year (evidence up to the year before)'); ax.set_ylabel('DisAtlas advantage over\nweighted Adamic–Adar (points)')
    ax.set_title(f"Across {NDRUG:,} drugs, DisAtlas is ahead in {NCI} of {NYRS} years (95% CI above zero)",loc='left',fontsize=7.5,fontweight='bold')
def draw_count_strip(ax, cen):
    ax.bar(cen.yr, cen.n_drugs, width=.7, color=CD, alpha=.35, lw=0)
    k=lambda n: f"{n/1000:.1f}k" if n>=1000 else f"{n}"
    for _,r in cen.iterrows():
        ax.text(r.yr, r.n_drugs+40, k(int(r.n_drugs)), ha='center', va='bottom', fontsize=4.3, color=CD)
        ax.text(r.yr, -120, k(int(r.n_events)), ha='center', va='top', fontsize=4.3, color=GR)
    ax.set_ylim(-520, 1750); ax.axis('off')
    ax.text(2021.9, 700, 'unique drugs', ha='left', va='center', fontsize=5, color=CD)
    ax.text(2021.9, -330, 'drug–disease tests', ha='left', va='center', fontsize=5, color=GR)
def build_figure(Ty, pooled, cen=None):
    fig=plt.figure(figsize=(7.2,6.0))
    gs=fig.add_gridspec(2,1,height_ratios=[3.6,2.2],hspace=.10)
    axa=fig.add_subplot(gs[0]); draw_panel_a(axa, CASE)
    sub=gs[1].subgridspec(2,3,width_ratios=[.07,.87,.06],height_ratios=[.42,1],hspace=.02)
    axs=fig.add_subplot(sub[0,1]); axb=fig.add_subplot(sub[1,1],sharex=axs)
    draw_panel_c(axb,Ty,pooled); t=axb.get_title(loc='left') or axb.get_title(); axb.set_title(''); axb.set_title('',loc='left')
    axb.spines[['top','right']].set_visible(False); axb.tick_params(labelsize=6); axb.xaxis.label.set_size(6.5); axb.yaxis.label.set_size(6.5)
    axb.set_xticks(range(1996,2022,2)); axb.set_xlim(1995.2,2024.5)
    draw_count_strip(axs,cen); axs.set_title(t,fontsize=8,fontweight='bold',loc='left',pad=4)
    fig.canvas.draw(); r_=fig.canvas.get_renderer(); inv=fig.transFigure.inverted()
    ya=inv.transform((0,axa.get_tightbbox(r_).y1))[1]; yb=inv.transform((0,axs.title.get_window_extent(r_).y1))[1]
    fig.text(0.07,ya,'a',fontsize=9,fontweight='bold',va='top',ha='left')
    fig.text(0.07,yb,'b',fontsize=9,fontweight='bold',va='top',ha='left')
    return fig


fig = build_figure(Ty, pooled, cen)
os.makedirs(OUT, exist_ok=True)
fig.savefig(f"{OUT}/fig7_prospective.png", dpi=int(os.environ.get("FIG_DPI", "300")), bbox_inches="tight", facecolor="white")
fig.savefig(f"{OUT}/fig7_prospective.pdf", bbox_inches="tight", facecolor="white")
from disatlas.parity import write_values
# stored deltas are weighted Adamic-Adar minus DisAtlas; values are written as the DisAtlas advantage (positive = better)
_lo, _hi = sorted([-pooled["lo"], -pooled["hi"]])
_V = {"n_drugs_pooled": NDRUG, "pooled_advantage": -pooled["d"], "pooled_advantage_lo": _lo, "pooled_advantage_hi": _hi}
for _r in Ty.itertuples(): _V[f"advantage|{_r.yr}"] = -_r.dE1
write_values(f"{OUT}/fig7_prospective_values.tsv", _V)
print("wrote", OUT)
