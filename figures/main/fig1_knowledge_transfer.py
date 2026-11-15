"""Main Figure 1 (rare and common diseases: genetic evidence and shared genes, Open Targets Platform 26.09) and
Supplementary Figure S1 (rare-disease statistics: a = prevalence panel, b = sources of genetic evidence).

Draws from the Figure 1 tables in data/derived/fig1_open_targets/ with disatlas.fig1_open_targets (subcommand
`figure`); the tables themselves are built from an Open Targets Platform release with the `tables` subcommand of the
same module (see its docstring).
Inputs (data/derived/fig1_open_targets/):
  fig1a_histogram_bins_OT26.09.csv   panel a histogram (% of diseases per 0.05 bin of the mean genetic-association score)
  fig1_summary_OT26.09.json          panel b gene counts and the group sizes
  si_source_shares_OT26.09.csv       Supplementary Figure S1b source shares
  examples_26.09.json                panel c/d example boxes (rare disease, shared genes, common disease)
  s1a_prevalence_panel.png           Supplementary Figure S1a, the prevalence panel used unchanged as an image
Writes figures/main/fig1_knowledge_transfer.{png,pdf} and figures/SI/rare_disease_stats.{png,pdf}.

    python figures/main/fig1_knowledge_transfer.py
"""
import os, sys, shutil, tempfile, types
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
import disatlas.deterministic  # noqa: F401  (byte-reproducible PDF output)
from disatlas.fig1_open_targets import cmd_figure

X = os.path.join(ROOT, "data", "derived", "fig1_open_targets")

if __name__ == "__main__":
    tmp = tempfile.mkdtemp(prefix="disatlas_fig1_")
    cmd_figure(types.SimpleNamespace(tables=X, release="26.09", examples=os.path.join(X, "examples_26.09.json"),
                                     s1a=os.path.join(X, "s1a_prevalence_panel.png"), out=tmp))
    for e in ("png", "pdf"):
        shutil.copyfile(os.path.join(tmp, f"knowledge_transfer.{e}"), os.path.join(HERE, f"fig1_knowledge_transfer.{e}"))
        shutil.copyfile(os.path.join(tmp, f"si_S1_composite.{e}"),
                        os.path.join(ROOT, "figures", "SI", f"rare_disease_stats.{e}"))
    shutil.rmtree(tmp)
    print("wrote", os.path.join(HERE, "fig1_knowledge_transfer.png"), os.path.join(ROOT, "figures", "SI", "rare_disease_stats.png"))
