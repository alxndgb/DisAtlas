"""Rebuild every main-text and Appendix table and check each cell against the manuscript.

Run from the repository root:  python tables/build_all.py
"""
import runpy
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SCRIPTS = ["main/table_protocols.py",
           "SI/table_SI_input_dis_gene_sources.py",
           "SI/table_SI_input_dis_symptom_sources.py",
           "SI/table_SI_ontology_baseline.py",
           "SI/table_SI_cluster_controls.py",
           "SI/table_SI_rare_common_strata.py",
           "SI/table_SI_sensitivity_post_test.py",
           "SI/table_SI_fig3a_trials.py"]

if __name__ == "__main__":
    for s in SCRIPTS:
        runpy.run_path(str(HERE / s), run_name="__main__")
