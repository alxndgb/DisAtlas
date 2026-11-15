"""Rebuild data/derived/degree_distributions/powerlaw_summary.tsv (Supplementary Figure SI_degree_distributions):
discrete power-law fit and goodness-of-fit test of the pooled degree distribution of each MechDist graph, by the
method of Clauset, Shalizi & Newman (2009, SIAM Rev 51:661) with the `powerlaw` package (Alstott et al. 2014).
Port of the authors' revision/response/_degree_distributions/powerlaw_csn.py. Not a model score: input is
data/derived/degree_distributions/degree_per_disease.tsv (deg_<graph>, diseases with degree >= 1 are fitted).

Needs the optional extra requirements-cluster.txt (powerlaw==2.0.0).

    python -m disatlas.rebuild.powerlaw_csn observed OUTDIR
    python -m disatlas.rebuild.powerlaw_csn boot_task OUTDIR TASK NCHUNK     (SLURM array, TASK in 0 .. 6*NCHUNK-1)
    python -m disatlas.rebuild.powerlaw_csn collect OUTDIR                   -> OUTDIR/powerlaw_summary.tsv

Fit settings: "free" = xmin by KS over all observed values; "floor" = xmin restricted above the construction floor
(broad 3, specific and combined 50). Alpha range widened to [1, 10] (powerlaw 2.0.0 defaults to [0, 3]).
Bootstrap (CSN 2009 sec. 4.1): 1,000 replicates per setting; replicate r of setting c uses
default_rng(20261101 + 100000 * index(c) + r); n_tail* ~ Binomial(n, n_tail / n) draws from the fitted discrete power
law (exact inverse CDF through the Hurwitz zeta function up to 1e6), the rest resampled from values below xmin;
xmin and alpha re-estimated on each replicate. p_GoF = share of replicates with KS D >= observed D.
Test: tests/test_rebuild_powerlaw_csn.py (observed fit always; full bootstrap from a sharded run, DISATLAS_POWERLAW_BOOT).
"""
import json, os, sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[3]
DD = ROOT / "data/derived/degree_distributions"
GRAPHS = ["broad", "specific", "combined"]
FLOOR = {"broad": 3, "specific": 50, "combined": 50}
CONFIGS = [(g, m) for g in GRAPHS for m in ("free", "floor")]
SEED0 = 20261101; NREP = 1000; KTAB = 10 ** 6
ALTS = ["lognormal", "exponential", "stretched_exponential", "truncated_power_law"]
PRANGE = {"alpha": [1.0, 10.0]}


def _pl():
    import powerlaw
    powerlaw.PARALLEL_ENABLE = False
    return powerlaw


def degrees(g):
    d = pd.read_csv(DD / "degree_per_disease.tsv", sep="\t")[f"deg_{g}"].to_numpy()
    return d, d[d >= 1].astype(float)


def fit(x, lo=None):
    pl = _pl()
    if lo is None:
        return pl.Fit(x, discrete=True, verbose=0, parameter_ranges=PRANGE)
    return pl.Fit(x, discrete=True, verbose=0, parameter_ranges=PRANGE, xmin=(float(lo), float(x.max())))


def lo_of(g, mode):
    return None if mode == "free" else FLOOR[g] + 1


def pl_sampler(alpha, xmin):
    from scipy.special import zeta
    k = np.arange(int(xmin), KTAB + 1, dtype=float)
    surv = zeta(alpha, k) / zeta(alpha, xmin); surv_desc = surv[::-1]
    def draw(rng, m):
        u = rng.random(m); out = np.empty(m); inside = u >= surv[-1]
        out[inside] = k[::-1][np.searchsorted(surv_desc, u[inside], side="left")]
        far = ~inside
        out[far] = np.floor((xmin - 0.5) * (1 - u[far]) ** (-1 / (alpha - 1)) + 0.5)
        return out
    return draw


def observed(out):
    res = {}
    for g, mode in CONFIGS:
        d, x = degrees(g); f = fit(x, lo_of(g, mode)); p = f.power_law
        r = dict(graph=g, fit=mode, xmin_lower_bound=lo_of(g, mode), n_diseases=int(len(d)), n_degree0=int((d == 0).sum()),
                 n_fit=int(len(x)), alpha=float(p.alpha), sigma=float(p.sigma), xmin=float(f.xmin), D=float(p.D),
                 n_tail=int((x >= f.xmin).sum()))
        for a in ALTS:
            R, pv = f.distribution_compare("power_law", a, normalized_ratio=True); r[f"LR_{a}_R"] = float(R); r[f"LR_{a}_p"] = float(pv)
        res[f"{g}:{mode}"] = r
    Path(out).mkdir(parents=True, exist_ok=True); json.dump(res, open(Path(out) / "powerlaw_observed.json", "w"), indent=1)
    return res


def boot(out, cfg, chunk, nchunk):
    g, mode = cfg.split(":"); lo = lo_of(g, mode)
    obs = json.load(open(Path(out) / "powerlaw_observed.json"))[cfg]
    _, x = degrees(g); n, xmin, alpha = len(x), obs["xmin"], obs["alpha"]
    below = x[x < xmin]; p_tail = obs["n_tail"] / n; draw = pl_sampler(alpha, xmin)
    (Path(out) / "boot").mkdir(parents=True, exist_ok=True); rows = []
    for r in [r for r in range(NREP) if r % nchunk == chunk]:
        seed = SEED0 + 100000 * CONFIGS.index((g, mode)) + r; rng = np.random.default_rng(seed)
        m = rng.binomial(n, p_tail)
        y = np.concatenate([draw(rng, m), rng.choice(below, n - m, replace=True) if len(below) else np.empty(0)])
        f = fit(y, lo)
        rows.append(f"{r}\t{seed}\t{m}\t{f.xmin:g}\t{f.power_law.alpha:.6f}\t{f.power_law.D:.8f}\n")
    p = Path(out) / "boot" / f"{g}_{mode}_{chunk:03d}.tsv"
    p.write_text("rep\tseed\tn_tail_star\txmin_star\talpha_star\tD_star\n" + "".join(rows))


def collect(out):
    obs = json.load(open(Path(out) / "powerlaw_observed.json")); rows = []
    for g, mode in CONFIGS:
        o = obs[f"{g}:{mode}"]; files = sorted((Path(out) / "boot").glob(f"{g}_{mode}_*.tsv"))
        b = pd.concat([pd.read_csv(p, sep="\t") for p in files]).drop_duplicates("rep")
        k = int((b.D_star >= o["D"]).sum()); nb = len(b); best = min(ALTS, key=lambda a: o[f"LR_{a}_R"])
        rows.append(dict(graph=g, fit=mode, xmin_lower_bound=o["xmin_lower_bound"], n_degree0=o["n_degree0"], n_diseases=o["n_diseases"],
                         n_fit=o["n_fit"], alpha=o["alpha"], alpha_mle=o["alpha"], sigma=o["sigma"], xmin=o["xmin"], n_tail=o["n_tail"], D=o["D"],
                         n_boot=nb, n_boot_D_ge_obs=k, p_gof=k / nb, p_gof_mc_se=float(np.sqrt(max(k, 0.5) / nb * (1 - k / nb) / nb)),
                         boot_xmin_median=float(b.xmin_star.median()), boot_alpha_median=float(b.alpha_star.median()),
                         best_rival=best, best_rival_R=o[f"LR_{best}_R"], best_rival_p=o[f"LR_{best}_p"],
                         **{f"LR_{a}_{s}": o[f"LR_{a}_{s}"] for a in ALTS for s in ("R", "p")}))
    S = pd.DataFrame(rows); S.to_csv(Path(out) / "powerlaw_summary.tsv", sep="\t", index=False, float_format="%.6g"); return S


if __name__ == "__main__":
    m, out = sys.argv[1], sys.argv[2]
    if m == "observed": observed(out)
    elif m == "boot_task":
        t, nch = int(sys.argv[3]), int(sys.argv[4]); g, mode = CONFIGS[t // nch]; boot(out, f"{g}:{mode}", t % nch, nch)
    elif m == "collect": print(collect(out).T.to_string())
