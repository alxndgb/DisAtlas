"""Numbers quoted in the SI cluster-validation, control-pair and rare-common sections.

Model-based quantities (control-pair link scores, percentile among all 84,766,710 disease pairs, best partner
rank, median link score, rheumatoid-arthritis case scores, Duchenne replacement check) are recomputed from the
released 5-seed model (data/zenodo/model). Everything else is read from the precomputed tables in
data/derived/si_cluster and data/derived/si_rare_common (see the README there). Writes SI_cluster_rarecommon.json.
"""
import json, gzip
from pathlib import Path
import numpy as np, pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from disatlas.sibio_common import derived, model_dir, write_numbers, SEEDS
from safetensors.numpy import load_file

N = {}
r3 = lambda x: round(float(x), 3)
# ---------------- model scoring ----------------
M = model_dir()
ids = [l.strip() for l in open(M / "nodes.txt") if l.strip()]
ix = {d: i for i, d in enumerate(ids)}
E, DEC = [], []
for s in SEEDS:
    z = np.load(M / "embeddings" / f"node2vec_s{s}.npz")
    assert list(z["ids"]) == ids
    E.append(z["X"].astype(np.float32))
    DEC.append(load_file(str(M / "decoders" / f"decoder_s{s}.safetensors")))

def score_rows(rows, cols=None):
    """ensemble probabilities for rows x cols (cols default: all diseases); returns (seed, rows, cols)"""
    out = []
    for X, d in zip(E, DEC):
        A = X[rows]; B = X if cols is None else X[cols]
        h = A[:, None, :] * B[None, :, :]
        h = np.maximum(h @ d["0.weight"].T + d["0.bias"], 0)
        z = (h @ d["3.weight"].T)[..., 0] + d["3.bias"][0]
        out.append(1 / (1 + np.exp(-z)))
    return np.stack(out)

n = len(ids); n_pairs = n * (n - 1) // 2
N["n_diseases"] = n; N["n_pairs_all"] = n_pairs
# all upper-triangle ensemble scores, in row blocks (histogram for percentiles + exact median)
allp = np.empty(n_pairs, np.float32); k = 0
for a in range(0, n, 64):
    rows = np.arange(a, min(a + 64, n))
    P = score_rows(rows).mean(0)
    for j, r in enumerate(rows):
        v = P[j, r + 1:]; allp[k:k + len(v)] = v; k += len(v)
assert k == n_pairs
allp.sort()
N["median_link_score_all_pairs"] = r3(np.median(allp))
pct = lambda p: 100.0 * np.searchsorted(allp, np.float32(p), side="left") / n_pairs

C = pd.read_csv(derived("si_cluster") / "SI_controls_table.tsv", sep="\t")
rows = []
for r in C.itertuples():
    i, j = ix[r.id_a], ix[r.id_b]
    Pi = score_rows(np.array([i, j])).mean(0)          # (2, n)
    p = float(Pi[0, j])
    def rank(row, self_i, other):
        v = np.delete(row, self_i); o = row[other]
        return int((v > o).sum() + 1)
    best = min(rank(Pi[0], i, j), rank(Pi[1], j, i))
    rows.append(dict(tag=r.tag, role=r.role, same=r.same_cluster_k55, score=round(p, 3), pct=round(pct(p), 1),
                     best_rank=best, table_score=r.disatlas_pair_score))
R = pd.DataFrame(rows)
N["controls"] = R[["tag", "role", "same", "score", "pct", "best_rank"]].to_dict("records")
pos, neg = R[R.role == "positive"], R[R.role == "negative"]
N["controls_n_pos"], N["controls_n_neg"] = len(pos), len(neg)
N["controls_pos_same_cluster"] = int((pos.same == "yes").sum())
N["controls_neg_separated"] = int((neg.same == "no").sum())
N["controls_split_pos_min_pct"] = float(pos[pos.same == "no"].pct.min())
N["controls_auroc"] = float(np.mean([[a > b for b in neg.score] for a in pos.score]))
N["controls_pos_top100"] = int((pos.best_rank <= 100).sum())
N["controls_pos_best_rank_range"] = [int(pos.best_rank.min()), int(pos.best_rank.max())]
N["controls_neg_best_rank_range"] = [int(neg.best_rank.min()), int(neg.best_rank.max())]
# Duchenne replacement check
DMD, DM1 = "MONDO:0010679", "MONDO:0008056"
Cand = pd.read_csv(derived("si_cluster") / "dmd_replacement_candidates.tsv", sep="\t")
el = Cand[Cand.eligible.astype(str) == "True"]
if DMD in ix:
    Pd = score_rows(np.array([ix[DMD]])).mean(0)[0]
    if DM1 in ix:
        N["dmd_dm1"] = [r3(Pd[ix[DM1]]), round(pct(Pd[ix[DM1]]), 1)]
    ep = [pct(Pd[ix[e]]) for e in el.id if e in ix]
    N["dmd_eligible"] = [len(ep), round(max(ep), 2)]
# Rheumatoid arthritis vs the IIM / overlap group (per-seed range)
# RA (EFO:0000685) and the IIM / overlap group, mapped to the MONDO ids of the 13,021-disease index via MONDO xrefs
RA = "MONDO:0008383"
RA_GROUP = [("polymyositis", "ORPHA:732", "MONDO:0019127"), ("antisynthetase", "ORPHA:81", "MONDO:0019344"),
            ("dermatomyositis_mondo", "MONDO:0016367", "MONDO:0016367"), ("mctd", "ORPHA:809", "MONDO:0005854"),
            ("classical_dm", "ORPHA:645613", "MONDO:0957456"), ("overlap_myositis", "ORPHA:206572", "MONDO:0016099")]
ra = {}
for lab, query, d in RA_GROUP:
    if d not in ix: ra[lab] = None; continue
    s = score_rows(np.array([ix[RA]]), np.array([ix[d]]))[:, 0, 0]
    ra[lab] = dict(query=query, mean=round(float(s.mean()), 3), lo=round(float(s.min()), 3), hi=round(float(s.max()), 3))
N["ra_case"] = ra; N["ra_id_found"] = RA in ix

# ---------------- precomputed tables ----------------
V = json.load(open(derived("si_cluster") / "sifu_cluster_13021_values.json"))
RL = json.load(open(derived("si_cluster") / "raw_vs_learned_silhouette_k55.json"))
N["sil_learned"], N["sil_cosine"], N["ari_cos_vs_learned"] = r3(RL["sil_learned"]), r3(RL["sil_raw_own_clusters"]), round(RL["ari_raw_vs_learned"], 2)
KB = pd.read_csv(derived("si_cluster") / "k_selection_boot.tsv.gz", sep="\t")
N["ksel_kinds"] = KB.kind.value_counts().to_dict()
KB = KB[KB.kind == "boot"]
kcol, scol = "K", "silhouette"
g = KB.groupby(kcol)[scol]; mean, sd = g.mean(), g.std()
kb = mean.idxmax(); thr = mean.max() - sd[kb]
N["ksel_best_mean_K"] = int(kb); N["ksel_best_mean"] = r3(mean.max()); N["ksel_best_sd"] = r3(sd[kb])
N["ksel_within_1sd"] = [int(x) for x in mean.index[mean >= thr]]
N["ksel_chosen_K"] = int(min(N["ksel_within_1sd"]))
q = KB[KB[kcol] == 55][scol]
N["ksel_K55"] = [r3(q.mean()), r3(q.quantile(.025)), r3(q.quantile(.975))]
N["ksel_K_range"] = [int(mean.index.min()), int(mean.index.max())]
N["ksel_n_subsamples"] = int(KB.groupby(kcol).size().max())
I = V["internal"]
N["internal"] = dict(sil=round(I["silhouette_all"], 2), within=round(I["median_within"], 2), nearest_other=round(I["median_nearest_other"], 2),
                     share_neg_sil=round(100 * I["share_diseases_negative_silhouette"], 1), n_clusters_neg=I["n_clusters_negative_mean_silhouette"],
                     size_range=[I["min_cluster_size"], I["max_cluster_size"]])
X = V["m3_external"]
N["external"] = dict(n_categories=X["n_categories"], n_with_category=X["diseases_with_category"],
                     same_cluster_share=round(100 * X["pair_share_within_cluster"]), random_share=round(100 * X["pair_share_random"]),
                     perm_p=round(X["perm_p_within"], 3), n_perm=X["n_null_perm_pairs"],
                     clusters_enriched=X["clusters_with_enriched_category"], enriched_null_max=X["clusters_with_enriched_category_null_max"],
                     n_shuffles=X["n_null_perm_enrichment"], ami=round(X["ami"], 2), ami_null_max=round(X["ami_null_max"], 3))
Mx = V["mixing"]
N["membership"] = dict(n_clusters=Mx["n_clusters"], both=Mx["both_rare_and_common"], common_ge10=Mx["common_share_ge_10pct"],
                       n_rare=Mx["n_rare"], n_common=Mx["n_common"], rare_share=round(100 * Mx["n_rare"] / (Mx["n_rare"] + Mx["n_common"])))
cc = V["cocluster"]["test"]
N["cocluster_test"] = dict(ratio=round(cc["ratio"], 2), pos=round(100 * cc["cocl_pos"], 1), neg=round(100 * cc["cocl_neg"], 1), n_pos=cc["n_pos"], n_neg=cc["n_neg"])
N["layout"] = dict(spearman=r3(V["layout"]["spearman"]), n_pairs=V["layout"]["n_pairs"])
# five display layouts: of the 9 x 8 = 72 positive-negative control comparisons, the number in which the positive
# pair lies nearer on the map (layout_comparison_summary.tsv; Isomap is the released map)
LC = pd.read_csv(derived("si_cluster") / "layout_comparison_summary.tsv", sep="\t").set_index("space")
N["layout_controls"] = {m: int(round(72 * LC.at[m, "pos_nearer"])) for m in ["isomap", "pca", "cmds", "umap", "tsne"]}
assert N["layout_controls"]["isomap"] == 72
# rare-common strata (IxIDN-Matched + CC+ supplementary)
S = pd.read_csv(derived("si_rare_common") / "ccplus_strata.tsv", sep="\t")
st = {}
for r in S.itertuples():
    key = f"{r.benchmark}|{r.scope}|{r.stratum}"
    st[key] = dict(n=int(r.n_pairs), pos=int(r.n_pos), neg=int(r.n_neg), auroc=r3(r.auroc), ci=[r3(r.ci_lo), r3(r.ci_hi)],
                   d_stack=r3(r.delta_vs_stack), d_ci=[r3(r.delta_vs_stack_lo), r3(r.delta_vs_stack_hi)], added=int(r.n_added_neg))
N["strata"] = st
F4 = pd.read_csv(derived("si_rare_common") / "four.tsv", sep="\t")
# anchored groups: mean per-anchor AUROC, anchor bootstrap B = 2000, seed 4707 (as si_rare_common_transfer.py)
anch = {}
for r in F4.itertuples():
    v = pd.read_csv(derived("si_rare_common") / f"{r.group}.tsv", sep="\t").auroc.values
    assert len(v) == r.n_anchors and abs(v.mean() - r.mean_auroc) < 6e-5, r.group
    bs = v[np.random.default_rng(4707).integers(0, len(v), (2000, len(v)))].mean(1)
    anch[r.group] = dict(auroc=r3(v.mean()), ci=[r3(np.percentile(bs, 2.5)), r3(np.percentile(bs, 97.5))],
                         anchors=int(r.n_anchors), pairs=int(r.n_pairs), neg=int(r.n_neg))
N["anchored"] = anch
EN = pd.read_csv(derived("si_rare_common") / "enrichment_other_class_frac10.tsv", sep="\t")
N["enrichment_k10"] = {f"{r.model}|{r.anchor_class}": [r3(r.mean_frac_other_in_10), round(r.p_lower, 3)] for r in EN.itertuples()}
CV = pd.read_csv(derived("si_rare_common") / "coverage_by_k_summary.tsv", sep="\t")
N["coverage_k10"] = {f"{r.model}|{r.anchor_class}": [round(100 * r.mean), round(r.p_lower, 3)] for r in CV[CV.k == 10].itertuples()}
RM = json.load(open(derived("si_rare_common") / "rc_coverage_meta.json"))
N["median_first_common_rank_for_rare"] = RM["median_first_common_rank_for_rare"]; N["coverage_nperm"] = RM["nperm"]
_obs = CV[(CV.k == 10) & (CV.model == "observed") & (CV.anchor_class == "rare")]["mean"].iloc[0]
N["coverage_k10_rare_none_pct"] = round(100 * (1 - _obs))   # rare diseases with no common disease among 10 nearest
# documented-partner recovery (68 test rare-common pairs): two-way disease-cluster bootstrap, B = 2000, seed 0;
# each draw resamples diseases with replacement and weights a pair by w[anchor] * w[partner]. The draws do not depend
# on the method columns, so adding TxGNN leaves the other methods' values unchanged. TxGNN column: published checkpoint
# scored as in the IxIDN-Matched comparison, pairs it cannot map ranked last (partner_percentile_per_pair_txgnn.tsv).
PP = pd.read_csv(derived("si_rare_common") / "partner_percentile_per_pair.tsv", sep="\t")
TX = pd.read_csv(derived("si_rare_common") / "partner_percentile_per_pair_txgnn.tsv", sep="\t")
PP = PP.merge(TX, on=["direction", "anchor", "partner"], how="left", validate="1:1")
assert PP.TxGNN.notna().all()
meth = ["DisAtlas", "cosine", "common_neighbours", "adamic_adar", "partner_degree", "TxGNN"]


def recovery(PP):
    rec = {}
    for dirn, gq in PP.groupby("direction", sort=False):
        rng = np.random.default_rng(0)
        dis = np.unique(np.r_[gq.anchor.values, gq.partner.values]); pos_ = {d: i for i, d in enumerate(dis)}
        a = gq.anchor.map(pos_).values; b = gq.partner.map(pos_).values; X = gq[meth].values
        est = {m: float(gq[m].mean()) for m in meth}; BS = []
        for _ in range(2000):
            w = rng.multinomial(len(dis), np.ones(len(dis)) / len(dis)).astype(float); ww = w[a] * w[b]
            BS.append((ww[:, None] * X).sum(0) / ww.sum())
        BS = np.array(BS); bs = {m: BS[:, k] for k, m in enumerate(meth)}
        out = {m: r3(est[m]) for m in meth}
        for m in meth[1:]:
            d = np.array(bs["DisAtlas"]) - np.array(bs[m]); out[f"d_{m}"] = [r3(est["DisAtlas"] - est[m]), r3(np.percentile(d, 2.5)), r3(np.percentile(d, 97.5))]
        out["n"] = len(gq); out["n_diseases"] = len(dis); rec[dirn] = out
    return rec


N["recovery"] = recovery(PP)
# sensitivity: only the pairs whose two diseases TxGNN can map (a different pair set, so the draws are redone on it)
N["recovery_txgnn_mapped"] = recovery(PP[PP.TxGNN_pair_mapped])
# IxIDN negative-pool funnel (benchmark construction; data/derived/si_rare_common/README.md)
FN = pd.read_csv(derived("si_rare_common") / "ixidn_negative_pool_funnel.tsv", sep="\t").set_index("stratum")
N["eligible_share_pct"] = {k: round(100 * r.negative_pool / r.all_pairs, 1) for k, r in FN.iterrows()}
N["negative_pool_funnel"] = {k: {c: int(v) for c, v in r.items()} for k, r in FN.iterrows()}
# headline test-split AUROC as quoted in the SI text (main-text disease-cluster bootstrap)
HT = pd.read_csv(derived("si_rare_common") / "headline_test_auroc.tsv", sep="\t").iloc[0]
N["headline_test"] = dict(auroc=r3(HT.auroc), ci=[r3(HT.ci_lo), r3(HT.ci_hi)], n=int(HT.n_pairs), pos=int(HT.n_pos))
N["n_scored_partners"] = int(N["n_diseases"]) - 1   # partners each disease ranks the other among
write_numbers(__file__, N, "SI_cluster_rarecommon")
print(len(N), "groups")
