"""Main Figure 4: how DisAtlas embedding cosines can be traced to shared evidence (gene-anchored model, cosine readout).

Panels:
  a test-split cosine distributions of IxIDN-Matched co-investigated and matched control pairs, tune-split threshold
  b test-split ROC of the cosine (AUROC with 95% disease-cluster bootstrap CI) and the operating point
  c cosine against the two-step walk mass two diseases share through neighbours (all IxIDN pairs, binned)
  d evidence trace of aHUS-PNH: top 8 evidence items, top 5 genes, top 3 shared neighbour diseases
  e share of walk mass by evidence type for the four case pairs
  f item-level traces of NMOSD-PNH, Crohn-UC and RA-polymyositis (top 5 items, 3 genes, 2 shared neighbours)
Displayed intermediates are ordered by share of walk mass x specificity weight ln(13,021 / degree); infectious diseases
and infectious-disease pathways follow the others (display only; no score changes; Supplementary Tables S22-S23).
Hub neighbour diseases (top 1% by degree) are not shown and Orphanet nodes identical to a listed MONDO node are merged.
Inputs: data/derived/fig4_evidence_trace/ (test_split_cosine.tsv, roc_cosine_test.tsv, m2_vs_cos_bins.tsv, fig4_values.json,
trace_displayed_items.tsv, trace_evidence_type_shares.tsv); only the displayed intermediates and per-type sums are shipped.

    python figures/main/fig4_evidence_trace.py
Writes fig4_evidence_trace.png/.pdf next to this script.
"""
import os, sys, json, textwrap, numpy as np, pandas as pd, matplotlib as mpl
mpl.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, os.path.join(REPO, 'src'))
import disatlas.deterministic  # noqa: F401  (date-free PDF metadata)
D = os.path.join(REPO, 'data', 'derived', 'fig4_evidence_trace')

# figure-style mechanics (size ladder 7 / 6 / 5.5, open frame, Type-42 fonts, 300 dpi)
S_BASE, S_SEC, S_TICK = 7, 6, 5.5
mpl.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'Nimbus Sans', 'DejaVu Sans'],
                     'font.size': S_BASE, 'axes.titlesize': S_BASE, 'axes.labelsize': S_BASE, 'axes.titlelocation': 'left',
                     'xtick.labelsize': S_TICK, 'ytick.labelsize': S_TICK, 'legend.fontsize': S_SEC, 'legend.frameon': False,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.6, 'lines.linewidth': 1.2,
                     'patch.linewidth': 0.6, 'xtick.major.size': 3, 'ytick.major.size': 3, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'pdf.fonttype': 42, 'ps.fonttype': 42, 'savefig.dpi': 300})
PLAB = {'aHUS-PNH': 'aHUS – PNH', 'NMOSD-PNH': 'NMOSD – PNH', 'Crohn-UC': 'Crohn – UC', 'RA-polymyositis': 'RA – polymyositis'}
C_BROAD = '#1b9e77'
C_POS, C_NEG = '#d95f02', '#7570b3'          # co-investigated vs matched control pairs
C_COS, C_GREY = '#111111', '#9a9a9a'

V = json.load(open(os.path.join(D, 'fig4_values.json')))
TE = pd.read_csv(os.path.join(D, 'test_split_cosine.tsv'), sep='\t')
ROC = pd.read_csv(os.path.join(D, 'roc_cosine_test.tsv'), sep='\t')
MB = pd.read_csv(os.path.join(D, 'm2_vs_cos_bins.tsv'), sep='\t')
TI = pd.read_csv(os.path.join(D, 'trace_displayed_items.tsv'), sep='\t')
GS = pd.read_csv(os.path.join(D, 'trace_evidence_type_shares.tsv'), sep='\t')
assert len(TE) == V['test_n_pos'] + V['test_n_neg'] and (TE.y == 1).sum() == V['test_n_pos']
SHOWN = TI[TI.shown_specificity].sort_values('rank_specificity', kind='stable')

CASE = ['aHUS-PNH', 'NMOSD-PNH', 'Crohn-UC', 'RA-polymyositis']
F_SHORT = {'Immunoregulatory interactions between a Lymphoid and a non-Lymphoid cell': 'Lymphoid–non-lymphoid interactions',
           'Translocation of ZAP-70 to Immunological synapse': 'ZAP-70 to immunological synapse',
           'Phosphorylation of CD3 and TCR zeta chains': 'CD3/TCR ζ phosphorylation',
           'Differentiation of Circulating Monocytes': 'Monocyte differentiation',
           'Regulation of Complement cascade': 'Complement cascade regulation',
           'Mendelian susceptibility to mycobacterial diseases due to complete IL12B deficiency': 'MSMD, IL12B deficiency',
           'Interleukin-10 signaling': 'IL-10 signalling', 'Interferon gamma signaling': 'IFN-γ signalling'}   # panel f display labels
D_SHORT = {'hemolytic uremic syndrome, atypical, susceptibility to, 1': 'aHUS susceptibility, type 1',
           'atypical hemolytic-uremic syndrome with B factor anomaly': 'aHUS with factor B anomaly'}   # panel d neighbour display labels
cap1 = lambda x: x if x.startswith('aHUS') else x[:1].upper() + x[1:]   # capitalise the first letter; aHUS keeps its lower-case a
DISP = []   # (pair, panel, list, rank, label, display label) of every displayed intermediate

def shown(pair, lst):
    return SHOWN[(SHOWN.pair == pair) & (SHOWN.list == lst)]

GROUP = {'Reactome pathway': 'Pathway', 'Literature pathway': 'Pathway', 'Literature GO process': 'GO process',
         'GO process (pathology)': 'GO process', 'Mechanism node (curated)': 'Mechanism node', 'Mechanism node (automated)': 'Mechanism node',
         'Literature mechanism': 'Literature mechanism', 'Named gene set': 'Named gene set', 'HPO phenotype': 'Phenotype, cell or tissue',
         'Cell type (pathology)': 'Phenotype, cell or tissue', 'Anatomy (pathology)': 'Phenotype, cell or tissue',
         'Shared rare gene': 'Shared rare gene', 'Broad similarity': 'Broad similarity', 'Class-factor residual': 'Class-factor residual', 'Other': 'Other', 'Gene': 'Gene'}
GCOL = {'Broad similarity': C_BROAD, 'Pathway': '#2b6cb0', 'GO process': '#8fb8de', 'Mechanism node': '#7b3294',
        'Literature mechanism': '#c2a5cf', 'Named gene set': '#e08214', 'Shared rare gene': '#35978f', 'Phenotype, cell or tissue': '#b2182b',
        'Class-factor residual': '#bdbdbd', 'Other': '#bdbdbd', 'Gene': '#252525', 'Neighbour disease': '#9a9a9a'}
GORDER = ['Broad similarity', 'Class-factor residual', 'Pathway', 'Mechanism node', 'Literature mechanism', 'Named gene set', 'Shared rare gene',
          'GO process', 'Phenotype, cell or tissue', 'Other']

def short(lbl, n=44):
    lbl = str(lbl)
    import re as _re; lbl = _re.sub(r' \(literature pattern v?\d+\)', ', literature', lbl)
    for a, b in [(' (complement generic)', ''), (' (complement C5 MAC)', ''), (' (complement C3 proximal)', ' proximal'), (' (autoantibody)', ''),
                 ('complement-fixing autoantibody (IgG1/IgG3 to cell-surface antigen; C1q->C3->C5->MAC; tissue C3/C9neo deposition; CDC)', 'complement-fixing autoantibody'),
                 ('broad channel similarity (GO, HPA, HPO, Pfam)', 'broad similarity (GO, HPA, HPO, Pfam)')]:
        lbl = lbl.replace(a, b)
    return lbl

# ---------------- panels ----------------
def p_hist(ax):
    t = V['threshold_youden_tune']; bins = np.linspace(-0.4, 1.0, 36)
    ax.set_ylim(0, 95)
    ax.hist(TE.c[TE.y == 0], bins=bins, color=C_NEG, alpha=0.6, label=f"Matched control (n = {V['test_n_neg']})")
    ax.hist(TE.c[TE.y == 1], bins=bins, color=C_POS, alpha=0.6, label=f"Co-investigated (n = {V['test_n_pos']})")
    ax.axvline(t, color=C_COS, lw=0.7, ls='--')
    ax.text(t + 0.025, 93, f"threshold {t:.2f}\nsensitivity {V['test_sensitivity']:.2f}\nspecificity {V['test_specificity']:.2f}",
            fontsize=S_SEC, va='top')
    ax.set_xlabel('Embedding cosine'); ax.set_ylabel('Test-split pairs')
    ax.legend(loc='upper right', handlelength=1, borderaxespad=0.2, bbox_to_anchor=(1.0, 0.70))
    ax.set_title('Co-investigated pairs have higher cosines than matched controls')

def p_roc(ax):
    ax.plot([0, 1], [0, 1], color=C_GREY, lw=0.5, ls=':')
    ax.plot(ROC.fpr, ROC.tpr, color=C_COS, lw=1.2,
            label=f"Cosine {V['cosine_test_auroc']:.3f} ({V['cosine_test_auroc_lo']:.3f}–{V['cosine_test_auroc_hi']:.3f})")
    ax.scatter([1 - V['test_specificity']], [V['test_sensitivity']], s=14, color=C_COS, zorder=4, label='Operating point')
    ax.set_xlim(-0.01, 1.01); ax.set_ylim(-0.01, 1.02); ax.set_aspect('equal')
    ax.set_xticks([0, 0.5, 1]); ax.set_yticks([0, 0.5, 1])
    ax.set_xlabel('False positive rate'); ax.set_ylabel('True positive rate')
    ax.legend(loc='lower right', fontsize=S_SEC, title='AUROC (95% CI)', title_fontsize=S_SEC, bbox_to_anchor=(1.04, 0),
              handlelength=1.2, borderaxespad=0.1, alignment='left')
    ax.set_title('The cosine separates test-split pairs')

def p_m2(ax):
    z = MB[MB.bin >= 0]
    ax.fill_between(z.log10M2_mid, z.cos_q25, z.cos_q75, color='#cccccc', lw=0)
    ax.plot(z.log10M2_mid, z.cos_med, color=C_COS, lw=1.0, marker='o', ms=2)
    z0 = MB[MB.bin < 0].iloc[0]; x0 = z.log10M2_mid.min() - 0.9
    ax.errorbar([x0], [z0.cos_med], yerr=[[z0.cos_med - z0.cos_q25], [z0.cos_q75 - z0.cos_med]], fmt='o', ms=2.5, color=C_GREY, lw=0.8)
    ax.text(x0 - 0.3, z0.cos_q75 + 0.05, 'no shared\nneighbour', ha='left', va='bottom', fontsize=S_TICK, color='#555')
    ax.set_xlabel('Shared walk mass (log10)'); ax.set_ylabel('Embedding cosine')
    ax.set_xlim(x0 - 0.6, z.log10M2_mid.max() + 0.3)
    ax.text(0.03, 0.97, f"Spearman {V['spearman_logM2_cos']:.2f}\n{V['n_m2_pos']:,} of {V['n_pairs']:,} pairs", transform=ax.transAxes, va='top', fontsize=S_SEC)
    ax.set_title('The cosine follows evidence\nshared through neighbours')

def p_trace(ax, pair='aHUS-PNH', k=8):
    mod, gen, nb = shown(pair, 'module'), shown(pair, 'gene'), shown(pair, 'neighbour')
    assert len(mod) == k and len(gen) == 5 and len(nb) == 3
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis('off')
    xa, xb, ya = 0.07, 0.93, 0.62
    ys = np.linspace(0.93, 0.33, k); wmax = mod.frac_of_M2.max()
    ax.text(0.5, ys[0] + 0.045, f'Top {len(mod)} items shown', ha='center', va='bottom', fontsize=S_TICK, color='#555')
    for i, (y, r) in enumerate(zip(ys, mod.itertuples()), 1):
        c = GCOL[GROUP[r.type]]; lw = 0.4 + 3.2 * r.frac_of_M2 / wmax
        for x0 in (xa, xb): ax.plot([x0 + (0.05 if x0 < 0.5 else -0.05), 0.5], [ya, y], color=c, lw=lw, alpha=0.75, solid_capstyle='butt', zorder=1)
        lb_ = short(r.label) + (', automated node' if r.type == 'Mechanism node (automated)' else '')
        DISP.append((pair, 'd', 'module', i, r.label, lb_))
        ax.text(0.5, y, f'{lb_} · {100 * r.frac_of_M2:.1f}%', ha='center', va='center', fontsize=S_TICK, zorder=3,
                bbox=dict(boxstyle='round,pad=0.22', fc='white', ec=c, lw=0.7))
    for x0, lb in [(xa, pair.split('-')[0]), (xb, pair.split('-')[1])]:
        ax.text(x0, ya, lb, ha='center', va='center', fontsize=S_BASE, fontweight='bold', zorder=3,
                bbox=dict(boxstyle='round,pad=0.35', fc='#f2f2f2', ec='#666666', lw=0.5))
    gs = [100 * r.frac_of_M2 for r in gen.itertuples()]; ns = [100 * r.frac_of_M2 for r in nb.itertuples()]
    gl = [r.label for r in gen.itertuples()]; nl = [cap1(D_SHORT.get(r.label, short(r.label))) for r in nb.itertuples()]
    DISP.extend((pair, 'd', 'gene', i, l, l) for i, l in enumerate(gl, 1)); DISP.extend((pair, 'd', 'neighbour', i, r.label, n) for i, (r, n) in enumerate(zip(nb.itertuples(), nl), 1))
    for yy, lead_, body, col in [(0.17, f'Top {len(gl)} genes shown', ' · '.join(f'{l} {s:.1f}%' for l, s in zip(gl, gs)), GCOL['Gene']),
                                 (0.08, f"Top {len(nl)} of {V['n_shared_neighbours'][pair]} shared neighbours", ' · '.join(f'{l} {s:.1f}%' for l, s in zip(nl, ns)), '#555')]:
        t0 = ax.text(0.0, yy, lead_, fontsize=S_SEC, va='center', color='#777')
        ax.annotate(body, xy=(1, 0.5), xycoords=t0, xytext=(8, 0), textcoords='offset points', va='center', ha='left', fontsize=S_SEC, color=col)
    ax.set_title('aHUS and PNH share complement evidence through neighbouring diseases')

def p_comp(ax):
    y = np.arange(len(CASE))[::-1]
    for yy, p in zip(y, CASE):
        v = GS[GS.pair == p].set_index('group').share; left = 0
        for gr in GORDER:
            if gr in v and v[gr] > 0:
                ax.barh(yy, v[gr], left=left, color=GCOL[gr], height=0.62, lw=0); left += v[gr]
    ax.set_yticks(y); ax.set_yticklabels([PLAB[p] for p in CASE]); ax.tick_params(axis='y', length=0)
    ax.set_xlim(0, 1); ax.set_xticks([0, 0.25, 0.5, 0.75, 1]); ax.set_xticklabels(['0', '25%', '50%', '75%', '100%']); ax.set_xlabel('Share of walk mass')
    ax.set_title('Evidence types behind\neach pair')

KF = {}
def p_items(ax, pair, km=5, kg=3, kn=2):
    mod, gen, nb = shown(pair, 'module'), shown(pair, 'gene'), shown(pair, 'neighbour')
    assert (len(mod), len(gen), len(nb)) == (km, kg, kn), pair
    fl = lambda r, suf='': cap1(F_SHORT.get(r.label, short(r.label)) + suf)
    rows = [(fl(r, ', node' if r.type == 'Mechanism node (automated)' else ''), r.frac_of_M2, GCOL[GROUP[r.type]]) for r in mod.itertuples()]
    rows += [(r.label, r.frac_of_M2, GCOL['Gene']) for r in gen.itertuples()]
    rows += [(fl(r), r.frac_of_M2, GCOL['Neighbour disease']) for r in nb.itertuples()]
    rows = [(textwrap.fill(l, 36), v, c_) for l, v, c_ in rows]
    src = [('module', r.label) for r in mod.itertuples()] + [('gene', r.label) for r in gen.itertuples()] + [('neighbour', r.label) for r in nb.itertuples()]
    rk_ = list(range(1, len(mod) + 1)) + list(range(1, len(gen) + 1)) + list(range(1, len(nb) + 1))
    DISP.extend((pair, 'f', li, k_, lab, l.replace('\n', ' ')) for (li, lab), k_, (l, _, _) in zip(src, rk_, rows))
    KF[pair] = (len(mod), len(gen), len(nb))
    gap = [0] * km + [0.6] * kg + [1.2] * kn
    y = -np.arange(len(rows)) - np.array(gap)
    ax.barh(y, [v for _, v, _ in rows], color=[c for *_, c in rows], height=0.7, lw=0)
    ax.set_yticks(y); ax.set_yticklabels([l for l, _, _ in rows], fontsize=S_TICK, linespacing=0.85); ax.tick_params(axis='y', length=0)
    for t, (_, _, c) in zip(ax.get_yticklabels(), rows): t.set_color('#222' if c != GCOL['Gene'] else GCOL['Gene'])
    rk = list(range(1, km + 1))[:len(mod)] + list(range(1, kg + 1))[:len(gen)] + list(range(1, kn + 1))[:len(nb)]
    ax.set_xlabel('Share of walk mass'); ax.set_xlim(0, 0.13); ax.set_xticks([0, 0.05, 0.10]); ax.set_xticklabels(['0', '5%', '10%'])
    ax.set_title(PLAB[pair], fontsize=S_SEC, fontweight='bold')
    ax._grp = [(y[:km], 'items'), (y[km:km + kg], 'genes'), (y[km + kg:], 'neighbours')]

def legend_types(fig, rect):
    lx = fig.add_axes(rect); lx.axis('off')
    grs = [gr for gr in GORDER if gr not in ('Other',) and gr != 'Class-factor residual']
    hs = [Patch(color=GCOL[gr], lw=0) for gr in grs] + [Patch(color=GCOL['Gene'], lw=0), Patch(color=GCOL['Neighbour disease'], lw=0)]
    lx.legend(hs, grs + ['Gene', 'Shared neighbour disease'], ncol=5, loc='center left', fontsize=S_SEC, handlelength=1.0, columnspacing=0.9,
              handletextpad=0.35, borderaxespad=0, bbox_to_anchor=(0, 0.5))

def letter(fig, axes_list, l, y=None):
    r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    bbs = [inv.transform(a.get_tightbbox(r)) for a in axes_list]
    x0 = min(b[0][0] for b in bbs); y1 = max(b[1][1] for b in bbs)
    fig.text(max(x0 - 0.012, 0.004), y if y is not None else y1 + 0.004, l, fontsize=9, fontweight='bold', va='bottom', ha='left')

# ---------------- layout ----------------
W_MM, H_MM = 180, 240
def f_layout(fig, gg, x_lo=0.012, x_hi=0.992, gap=0.025, pad=0.008):
    # one set of 'Top k' brackets left of the first column; columns placed so no label touches the previous column or the brackets
    assert len(set(KF.values())) == 1, f'k differs between panel f columns: {KF}'
    km, kg, kn = next(iter(KF.values())); rr = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    ax0 = gg[0]; brk = []
    for (yy, nm), k_ in zip(ax0._grp, (km, kg, kn)):
        brk.append(fig.text(0, 0, f'Top {k_}\n{nm}', fontsize=S_TICK, color='#555', ha='right', va='center', linespacing=0.95))
    fig.canvas.draw(); wb = max(inv.transform(t.get_window_extent(rr))[1][0] - inv.transform(t.get_window_extent(rr))[0][0] for t in brk)
    ext = []
    for ax in gg:
        pos = ax.get_position(); tb = inv.transform(ax.get_tightbbox(rr))
        ext.append((pos.x0 - tb[0][0], tb[1][0] - pos.x1))
    zone = wb + 3 * pad + 0.006
    w = (x_hi - x_lo - zone - sum(l + r for l, r in ext) - gap * (len(gg) - 1)) / len(gg); assert w > 0.06, w
    x = x_lo + zone
    for ax, (l, r) in zip(gg, ext):
        pos = ax.get_position(); ax.set_position([x + l, pos.y0, w, pos.height]); x = x + l + w + r + gap
    fig.canvas.draw()
    ext2 = []
    for ax in gg:
        tb = inv.transform(ax.get_tightbbox(rr)); ext2.append((tb[0][0], tb[1][0]))
    gaps = [ext2[i + 1][0] - ext2[i][1] for i in range(len(gg) - 1)]; assert min(gaps) > 0.01, gaps
    xl = ext2[0][0] - pad   # bracket line x, left of the first column's labels
    for (yy, nm), t in zip(ax0._grp, brk):
        y0, y1 = (fig.transFigure.inverted().transform(ax0.transData.transform((0, v)))[1] for v in (yy.max() + 0.4, yy.min() - 0.4))
        fig.add_artist(Line2D([xl + 0.005, xl, xl, xl + 0.005], [y0, y0, y1, y1], transform=fig.transFigure, color='#555', lw=0.5))
        t.set_position((xl - pad, (y0 + y1) / 2))
    fig.canvas.draw()

def build(out):
    fig = plt.figure(figsize=(W_MM / 25.4, H_MM / 25.4))
    a = fig.add_axes([0.07, 0.80, 0.36, 0.15]); p_hist(a)
    b = fig.add_axes([0.50, 0.80, 0.19, 0.15]); p_roc(b)
    c = fig.add_axes([0.80, 0.80, 0.17, 0.125]); p_m2(c)
    d = fig.add_axes([0.0, 0.43, 0.62, 0.28]); p_trace(d)
    e = fig.add_axes([0.79, 0.475, 0.175, 0.245]); p_comp(e)
    legend_types(fig, [0.01, 0.392, 0.98, 0.022])
    gg = [fig.add_axes([x, 0.055, 0.115, 0.28]) for x in (0.215, 0.535, 0.865)]
    for ax, p in zip(gg, CASE[1:]): p_items(ax, p)
    fig.text(0.03, 0.362, 'Items behind the other case pairs', fontsize=S_BASE, va='bottom')
    fig.canvas.draw()
    f_layout(fig, gg)
    rows = [[[a], [b], [c]], [[d], [e]], [gg]]
    rr = fig.canvas.get_renderer(); inv = fig.transFigure.inverted(); L_ = iter('abcdefgh')
    for row in rows:
        ytop = max(inv.transform(x.get_tightbbox(rr))[1][1] for grp in row for x in grp)
        for grp in row:
            l = next(L_)
            if grp is gg: fig.text(0.004, 0.362, l, fontsize=9, fontweight='bold', va='bottom')
            elif grp[0] is d: fig.text(0.004, ytop + 0.004, l, fontsize=9, fontweight='bold', va='bottom')
            else: letter(fig, grp, l, y=ytop + 0.004)
    for ext in ['png', 'pdf']: fig.savefig(f'{out}.{ext}', dpi=300)
    plt.close(fig)

if __name__ == '__main__':
    build(os.path.join(HERE, 'fig4_evidence_trace'))
    # the display labels drawn must equal the shipped display_label column (Supplementary Table S18)
    G = pd.DataFrame(DISP, columns=['pair', 'panel', 'list', 'rank', 'label', 'display_label']).drop_duplicates()
    S = SHOWN.rename(columns={'rank_specificity': 'rank'})[['pair', 'panel', 'list', 'rank', 'label', 'display_label']]
    m = G.merge(S, on=['pair', 'panel', 'list', 'rank'], suffixes=('', '_shipped'))
    assert len(m) == len(G) == len(S) and (m.label == m.label_shipped).all() and (m.display_label == m.display_label_shipped).all()
    print('fig4_evidence_trace: displayed intermediates', len(G))
