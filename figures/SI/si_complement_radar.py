"""Panel b of the complement-milestones SI figure: radar of the DisAtlas cosine score (left) and its percentile among all
84,766,710 disease pairs (right; radial axis from the 90th percentile) for the six pairs of eculizumab's four approved
indications PNH (FDA 2007), aHUS (2011), gMG (2017) and NMOSD (2019); one spoke per pair.
Scores come only from the shipped five-seed mean cosine matrix (data/zenodo/model_outputs/cosine_matrix, seeds 4707-4711);
percentile = 100 * #(upper-triangle pairs with cosine <= c) / 84,766,710 (np.searchsorted(sorted, c, 'right')), the
all-pairs percentile used in figures/main/fig4_numbers.py. Imported by si_complement_milestones.py."""
import json, numpy as np, pandas as pd, matplotlib as mpl
from disatlas import paths

DIS = [("PNH", 2007), ("aHUS", 2011), ("gMG", 2017), ("NMOSD", 2019)]
MONDO = dict(PNH="MONDO:0100244", aHUS="MONDO:0016244", gMG="MONDO:0009688", NMOSD="MONDO:0019100")
EDGE = dict(PNH="#08306b", aHUS="#2171b5", gMG="#4292c6", NMOSD="#6baed6")
SCALE = dict(cos=("cosine", 0.0, 1.0, [0.25, 0.5, 0.75, 1.0], "cosine score"),
             pct=("percentile", 90.0, 100.0, [92.5, 95, 97.5, 100], "percentile of all disease pairs"))
N_PAIRS = 84_766_710


def pairs_table():
    """Six eculizumab-indication pairs: five-seed mean cosine and all-pairs percentile from the shipped matrix."""
    MO = paths.ZENODO / "model_outputs" / "cosine_matrix"
    ids = json.load(open(MO / "ids_13021.json")); ix = {d: i for i, d in enumerate(ids)}
    C = np.load(MO / "cosine_ens5_13021.npy", mmap_mode="r"); n = len(ids)
    ref = np.sort(np.concatenate([np.asarray(C[i, i + 1:], dtype=np.float32) for i in range(n - 1)])); assert ref.size == N_PAIRS, ref.size
    rows = []
    for i in range(4):
        for j in range(i + 1, 4):
            a, b = DIS[i][0], DIS[j][0]; c = np.float32(C[ix[MONDO[a]], ix[MONDO[b]]])
            rows.append(dict(pair=f"{a}-{b}", disease_a=a, id_a=MONDO[a], approval_a=DIS[i][1], disease_b=b, id_b=MONDO[b],
                             approval_b=DIS[j][1], cosine=round(float(c), 4),
                             percentile=round(100.0 * np.searchsorted(ref, c, "right") / N_PAIRS, 4)))
    return pd.DataFrame(rows)


def load(P):
    assert len(P) == 6, len(P)
    v = {}
    for r in P.itertuples():
        v[(r.disease_a, r.disease_b)] = v[(r.disease_b, r.disease_a)] = r
    return v


def _frame(ax, lo, hi, ticks, fmt):
    ax.set_ylim(lo, hi); ax.set_yticks(ticks); ax.set_yticklabels([fmt(t) for t in ticks], fontsize=5.5, color="0.35")
    ax.set_rlabel_position(210 if len(ax.get_xticks()) == 6 else 45)  # between the two lower-left spokes, clear of the value labels; ax.grid(color="0.85", lw=0.5); ax.spines["polar"].set_color("0.6"); ax.spines["polar"].set_linewidth(0.5)


def draw(ax, v, variant, scale, title=None):
    col, lo, hi, ticks, rlab = SCALE[scale]
    fmt = (lambda t: f"{t:g}") if scale == "cos" else (lambda t: f"{t:g}")
    val = lambda a, b: getattr(v[(a, b)], col)
    assert variant == "b", variant  # six pair spokes
    pairs = [(DIS[i], DIS[j]) for i in range(4) for j in range(i + 1, 4)]
    th = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    r = [val(a[0], b[0]) for a, b in pairs]
    ax.fill(np.r_[th, th[0]], np.r_[r, r[0]], color="#2166ac", alpha=0.15, zorder=2)
    ax.plot(np.r_[th, th[0]], np.r_[r, r[0]], color="#2166ac", lw=1.1, zorder=3)
    ax.plot(th, r, "o", ms=2.8, color="#2166ac", zorder=4)
    for t_, r_ in zip(th, r):  # value label on the radial side with more room: outward near the centre (clear of the
        phi = np.pi / 2 - t_     # radial tick labels), inward near the rim (clear of the spoke labels)
        sgn = 1 if (r_ - lo) / (hi - lo) < 0.6 else -1; dx, dy = sgn * 9 * np.cos(phi), sgn * 9 * np.sin(phi)
        ax.annotate(f"{r_:.3f}" if scale == "cos" else f"{r_:.1f}", (t_, r_), xytext=(dx, dy), textcoords="offset points",
                    ha="center", va="center", fontsize=5.5, color="#08306b", zorder=5,
                    bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none", alpha=0.8))
    ax.set_xticks(th); ax.set_xticklabels([f"{a[0]}–{b[0]}\n{a[1]}/{b[1]}" for a, b in pairs], fontsize=6.5)
    ax.set_theta_offset(np.pi / 2); ax.set_theta_direction(-1); ax.set_xlim(0, 2 * np.pi)
    ax.tick_params(axis="x", pad=1)
    _frame(ax, lo, hi, ticks, fmt)
    ax.set_title(title or rlab, fontsize=7, loc="center", pad=12)
