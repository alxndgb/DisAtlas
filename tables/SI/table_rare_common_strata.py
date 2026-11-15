"""Supplementary Table SI_rare_common_strata: AUROC of DisAtlas on IxIDN-Matched by disease-rarity pairing.

Rows: rare-rare, rare-common and common-common pairs of IxIDN-Matched, and the common-common stratum of the
supplementary benchmark with added matched negative pairs. AUROC, disease-cluster bootstrap 95% CI and the paired
difference to the combined matched-shortcut baselines (gene-anchored model, embedding-cosine score) are read from
data/derived/si_rare_common/ccplus_strata.tsv (the same disease-cluster bootstrap draws as Figure 5c).
Outputs next to this script: table_rare_common_strata.{md,tsv}
"""
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
T = pd.read_csv(ROOT / "data/derived/si_rare_common/ccplus_strata.tsv", sep="\t")
S = {f"{r.benchmark}|{r.scope}|{r.stratum}": dict(n=int(r.n_pairs), pos=int(r.n_pos), neg=int(r.n_neg), added=int(r.n_added_neg),
     auroc=r.auroc, ci=(r.ci_lo, r.ci_hi), d_stack=r.delta_vs_stack, d_ci=(r.delta_vs_stack_lo, r.delta_vs_stack_hi))
     for r in T.itertuples()}
f3 = lambda x: f"{x:.3f}"
sg = lambda x: ("+" if x >= 0 else "−") + f"{abs(x):.3f}"


def lohi(ci):
    return f"{f3(ci[0])}–{f3(ci[1])}"


def delta(s, extra=""):
    lo, hi = s["d_ci"]
    sig = lo > 0 or hi < 0
    return f"{sg(s['d_stack'])} ({sg(lo)} to {sg(hi)}{extra}{'' if sig else '; not significant'})"


m = {k.split("|")[2]: v for k, v in S.items() if k.startswith("IxIDN-Matched|All pairs|")}
cp = S["IxIDN-Matched-CC+|All pairs|common-common"]
cc = m["common-common"]
rows = [["Rare–rare", f"{m['rare-rare']['n']:,}", f3(m["rare-rare"]["auroc"]), lohi(m["rare-rare"]["ci"]), delta(m["rare-rare"])],
        ["Rare–common", f"{m['rare-common']['n']:,}", f3(m["rare-common"]["auroc"]), lohi(m["rare-common"]["ci"]), delta(m["rare-common"])],
        ["Common–common (IxIDN-Matched)", f"{cc['n']:,} ({cc['pos']:,} / {cc['neg']:,})", f3(cc["auroc"]), lohi(cc["ci"]),
         delta(cc, f"; {cc['neg']} negatives")],
        ["Common–common (supplementary, matched negatives)", f"{cp['added']} added negatives", f3(cp["auroc"]), lohi(cp["ci"]),
         "not significant" if cp["d_ci"][0] <= 0 <= cp["d_ci"][1] else delta(cp)]]
cols = ["Stratum", "Pairs", "AUROC", "95% CI", "DisAtlas − combined matched-shortcut baselines (95% CI)"]
pd.DataFrame(rows, columns=cols).to_csv(HERE / "table_rare_common_strata.tsv", sep="\t", index=False)
md = ["| " + " | ".join(cols) + " |", "|:--|:--|--:|--:|:--|"] + ["| " + " | ".join(r) + " |" for r in rows]
(HERE / "table_rare_common_strata.md").write_text("\n".join(md) + "\n", encoding="utf-8")
print(rows)
