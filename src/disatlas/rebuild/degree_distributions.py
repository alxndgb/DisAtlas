"""Rebuild data/derived/degree_distributions/{degree_summary.tsv, powerlaw_summary.tsv} (Supplementary Figure S10:
number of linked diseases (degree) in each graph of the MechDist network, with the power-law tail fits).

    python -m disatlas.rebuild.degree_distributions summary OUTDIR
    python -m disatlas.rebuild.degree_distributions observed OUTDIR
    python -m disatlas.rebuild.degree_distributions boot OUTDIR CONFIG CHUNK NCHUNK      # CONFIG = graph:fit, e.g. broad:free
    python -m disatlas.rebuild.degree_distributions collect OUTDIR

Input (earliest releasable input): data/derived/degree_distributions/degree_per_disease.tsv, the degree of each of the
13,021 diseases in the broad channel, the specific channel and the combined graph, with the rare / common class. The
498,155-edge combined graph the degrees were counted on is not shipped (the released graph of record is the later P_full
graph), so the degrees themselves are the input; the checks that tie them to the stated edge counts are in
tests/test_rebuild_degree_distributions.py.

summary  = degree_summary.tsv (per graph and group: isolated diseases, mean / median / quantiles, top-1% endpoint share,
           Mann-Whitney rare vs common), float format %.6g as stored.
observed = discrete power-law fit of the pooled degrees >= 1 by the method of Clauset, Shalizi & Newman (2009, SIAM Rev
           51:661) with the `powerlaw` package (Alstott et al. 2014; powerlaw 2.0.0, optional extra
           requirements-powerlaw.txt): alpha (MLE), sigma, xmin (KS-minimising), KS D, n_tail, and the normalised
           likelihood-ratio tests against lognormal, exponential, stretched-exponential and truncated power law.
           Two fit settings per graph: "free" (xmin chosen over all observed values) and "floor" (xmin strictly above the
           construction floor of the graph: 3 broad, 50 specific and combined).
boot     = semi-parametric bootstrap of the goodness of fit (CSN 2009 section 4.1): per replicate n_tail* ~ Binomial(n,
           n_tail / n) values from the fitted discrete power law (exact inverse CDF via the Hurwitz zeta function) and
           n - n_tail* observed values below xmin, refitted exactly as the data. Seed of replicate r of setting c:
           SEED0 + 100000 * index(c) + r. A replicate that exceeds FIT_TIMEOUT_S is written with status=timeout.
collect  = p_GoF = share of replicates with KS D >= observed D, written as powerlaw_summary.tsv and the per-replicate table.
"""
import json
import os
import signal
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import zeta
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[3]
DEG = ROOT / "data/derived/degree_distributions/degree_per_disease.tsv"
GRAPHS = ["broad", "specific", "combined"]
FLOOR = {"broad": 3, "specific": 50, "combined": 50}
CONFIGS = [(g, m) for g in GRAPHS for m in ("free", "floor")]
SEED0 = 20261101
NREP = 1000
ALTS = ["lognormal", "exponential", "stretched_exponential", "truncated_power_law"]
KTAB = 10 ** 6
PRANGE = {"alpha": [1.0, 10.0]}
FIT_TIMEOUT_S = 900


def load():
    return pd.read_csv(DEG, sep="\t")


def summary(P):
    cls = P.disease_class.to_numpy(); n = len(P)
    rare, common = cls == "rare", cls == "common"
    assert rare.sum() == 10489 and common.sum() == 2532 and (rare | common).all()
    groups = {"all": np.ones(n, bool), "rare": rare, "common": common}
    rows = []
    for g in GRAPHS:
        d = P[f"deg_{g}"].to_numpy(); L = int(d.sum() // 2); srt = np.sort(d)[::-1]
        for grp, m in groups.items():
            x = d[m]; nz = x[x > 0]
            r = dict(graph=g, group=grp, n_diseases=int(m.sum()), n_edges_graph=L, n_isolated=int((x == 0).sum()),
                     pct_isolated=100 * (x == 0).mean(), mean_degree=x.mean(), median_degree=float(np.median(x)),
                     q25=float(np.percentile(x, 25)), q75=float(np.percentile(x, 75)), p90=float(np.percentile(x, 90)),
                     p99=float(np.percentile(x, 99)), max_degree=int(x.max()),
                     median_degree_connected=float(np.median(nz)) if len(nz) else np.nan)
            if grp == "all":
                k1 = max(1, int(round(0.01 * n)))
                r["pct_endpoints_top1pct"] = 100 * srt[:k1].sum() / srt.sum()
                r["pct_endpoints_top10pct"] = 100 * srt[:int(round(0.1 * n))].sum() / srt.sum()
                r["n_top1pct"] = k1
            rows.append(r)
        mw = mannwhitneyu(d[rare], d[common], alternative="two-sided")
        rows.append(dict(graph=g, group="rare_vs_common", n_diseases=n, n_edges_graph=L,
                         auc_rare_gt_common=mw.statistic / (rare.sum() * common.sum()), mannwhitney_p=mw.pvalue,
                         median_ratio_common_over_rare=float(np.median(d[common])) / max(float(np.median(d[rare])), 1e-9)))
    return pd.DataFrame(rows)


def degrees(g):
    d = load()[f"deg_{g}"].to_numpy()
    return d, d[d >= 1].astype(float)


def fit(x, lo=None):
    import powerlaw
    powerlaw.PARALLEL_ENABLE = False
    if lo is None:
        return powerlaw.Fit(x, discrete=True, verbose=0, parameter_ranges=PRANGE)
    return powerlaw.Fit(x, discrete=True, verbose=0, parameter_ranges=PRANGE, xmin=(float(lo), float(x.max())))


def lo_of(g, mode):
    return None if mode == "free" else FLOOR[g] + 1


def pl_sampler(alpha, xmin):
    k = np.arange(int(xmin), KTAB + 1, dtype=float)
    surv = zeta(alpha, k) / zeta(alpha, xmin)
    surv_desc = surv[::-1]

    def draw(rng, m):
        u = rng.random(m)
        out = np.empty(m)
        inside = u >= surv[-1]
        idx = np.searchsorted(surv_desc, u[inside], side="left")
        out[inside] = k[::-1][idx]
        far = ~inside
        out[far] = np.floor((xmin - 0.5) * (1 - u[far]) ** (-1 / (alpha - 1)) + 0.5)
        return out
    return draw


def observed(configs=CONFIGS):
    res = {}
    for g, mode in configs:
        d, x = degrees(g)
        f = fit(x, lo_of(g, mode)); pl = f.power_law
        r = dict(graph=g, fit=mode, xmin_lower_bound=lo_of(g, mode), floor=FLOOR[g], n_diseases=int(len(d)),
                 n_degree0=int((d == 0).sum()), n_fit=int(len(x)), alpha=float(pl.alpha), sigma=float(pl.sigma),
                 xmin=float(f.xmin), D=float(pl.D), n_tail=int((x >= f.xmin).sum()), max_degree=int(x.max()))
        for a in ALTS:
            R, p = f.distribution_compare("power_law", a, normalized_ratio=True)
            r[f"LR_{a}_R"] = float(R); r[f"LR_{a}_p"] = float(p)
        res[f"{g}:{mode}"] = r
    return res


class _Timeout(Exception):
    pass


def _alarm(signum, frame):
    raise _Timeout()


def boot_replicate(cfg, obs, r):
    """One bootstrap replicate r of setting cfg ('graph:fit') given the observed fit dict."""
    g, mode = cfg.split(":"); lo = lo_of(g, mode)
    _, x = degrees(g)
    n, xmin, alpha = len(x), obs["xmin"], obs["alpha"]
    below = x[x < xmin]; p_tail = obs["n_tail"] / n
    seed = SEED0 + 100000 * CONFIGS.index((g, mode)) + r
    rng = np.random.default_rng(seed)
    m = rng.binomial(n, p_tail)
    y = np.concatenate([pl_sampler(alpha, xmin)(rng, m), rng.choice(below, n - m, replace=True) if len(below) else np.empty(0)])
    t0 = time.time(); signal.signal(signal.SIGALRM, _alarm); signal.alarm(FIT_TIMEOUT_S)
    try:
        f = fit(y, lo); signal.alarm(0)
        return dict(cfg=cfg, rep=r, seed=seed, n_tail_star=int(m), xmin_star=float(f.xmin), alpha_star=round(float(f.power_law.alpha), 6),
                    D_star=round(float(f.power_law.D), 8), status="ok", seconds=round(time.time() - t0, 1))
    except _Timeout:
        return dict(cfg=cfg, rep=r, seed=seed, n_tail_star=int(m), xmin_star=np.nan, alpha_star=np.nan, D_star=np.nan,
                    status="timeout", seconds=round(time.time() - t0, 1))


def boot(out, cfg, chunk, nchunk):
    obs = json.load(open(Path(out) / "powerlaw_observed.json"))[cfg]
    rows = []
    for r in range(NREP):
        if r % nchunk != chunk:
            continue
        rows.append(boot_replicate(cfg, obs, r)); print("rep", r, rows[-1]["status"], rows[-1]["seconds"], flush=True)
    (Path(out) / "boot").mkdir(exist_ok=True)
    pd.DataFrame(rows).to_csv(Path(out) / "boot" / f"{cfg.replace(':', '_')}_{chunk:03d}.tsv", sep="\t", index=False)


def collect_table(obs, B):
    """B: per-replicate table (cfg, rep, xmin_star, alpha_star, D_star, status). Returns the powerlaw_summary frame."""
    rows = []
    for g, mode in CONFIGS:
        o = obs[f"{g}:{mode}"]; b = B[B.cfg == f"{g}:{mode}"].drop_duplicates("rep")
        if not len(b):
            continue
        assert (b.status == "ok").all(), "timed-out replicates present"
        k = int((b.D_star >= o["D"]).sum()); nb = len(b)
        best = min(ALTS, key=lambda a: o[f"LR_{a}_R"])
        rows.append(dict(graph=g, fit=mode, xmin_lower_bound=o["xmin_lower_bound"], n_degree0=o["n_degree0"], n_diseases=o["n_diseases"],
                         n_fit=o["n_fit"], alpha=o["alpha"], alpha_mle=o["alpha"], sigma=o["sigma"], xmin=o["xmin"], n_tail=o["n_tail"], D=o["D"],
                         n_boot=nb, n_boot_D_ge_obs=k, p_gof=k / nb,
                         p_gof_mc_se=float(np.sqrt(max(k, 0.5) / nb * (1 - k / nb) / nb)) if nb else np.nan,
                         boot_xmin_median=float(b.xmin_star.median()), boot_alpha_median=float(b.alpha_star.median()),
                         best_rival=best, best_rival_R=o[f"LR_{best}_R"], best_rival_p=o[f"LR_{best}_p"],
                         **{f"LR_{a}_{s}": o[f"LR_{a}_{s}"] for a in ALTS for s in ("R", "p")}))
    return pd.DataFrame(rows)


def collect(out):
    out = Path(out)
    obs = json.load(open(out / "powerlaw_observed.json"))
    B = pd.concat([pd.read_csv(p, sep="\t") for p in sorted((out / "boot").glob("*.tsv"))])
    B.sort_values(["cfg", "rep"]).to_csv(out / "powerlaw_boot_replicates.tsv", sep="\t", index=False)
    S = collect_table(obs, B); S.to_csv(out / "powerlaw_summary.tsv", sep="\t", index=False, float_format="%.6g")
    print(S.T.to_string())


if __name__ == "__main__":
    m, out = sys.argv[1], Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
    if m == "summary":
        summary(load()).to_csv(out / "degree_summary.tsv", sep="\t", index=False, float_format="%.6g")
    elif m == "observed":
        json.dump(observed(), open(out / "powerlaw_observed.json", "w"), indent=1)
    elif m == "boot":
        boot(out, sys.argv[3], int(sys.argv[4]), int(sys.argv[5]))
    elif m == "boot_task":              # array index -> (config, chunk)
        t, nch = int(sys.argv[3]), int(sys.argv[4]); g, mode = CONFIGS[t // nch]; boot(out, f"{g}:{mode}", t % nch, nch)
    elif m == "collect":
        collect(out)
