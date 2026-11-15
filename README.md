# DisAtlas

Code and data to reproduce the results of the manuscript

> **A Mechanistic Disease Representation Framework for Therapeutic Knowledge Transfer Across Rare and Common Conditions**

and its Supplementary Information, including the trial co-investigation benchmarks **IxIDN+**, **IxIDN-** and
**IxIDN-Matched**. The scope is the revised manuscript and its Supplementary Information: every main and
Supplementary figure and table that is computed from data, and every number in the main text, is regenerated
by one command and checked by the test suite.

DisAtlas places 13,021 diseases (10,489 rare, 2,532 common) in one embedding space learned from a
disease-disease network (MechDist) with two evidence channels: a **broad** channel (shared gene function, phenotype,
cell-type expression and protein domains) and a **specific** channel (shared mechanism modules). A five-seed
node2vec ensemble embeds every disease; a pair is scored by the cosine similarity of its two disease vectors.

## Contents

1. [Manuscript to repository map](#manuscript-to-repository-map)
2. [Repository structure](#repository-structure)
3. [Quick start](#quick-start)
4. [Run modes](#run-modes)
5. [Data and licensing](#data-and-licensing)
6. [Tests and parity](#tests-and-parity)
7. [Troubleshooting](#troubleshooting)
8. [How to cite](#how-to-cite)
9. [Licence](#licence)

## Manuscript to repository map

Numbering follows the final revised manuscript and SI. Every script runs from the repository root without arguments, except Figures S20 and S21, which are drawn by the Figure 7 script with `--si hits10` / `--si weighted` (`./run_pipeline.sh` calls them in order). Each script writes its outputs next to itself; most also write a `*_values.tsv` / `*_values.json` / `*_numbers.json` file with every plotted or quoted number, which the tests read. All DisAtlas pair scores are embedding cosines (mean of the five per-seed cosines).

Links: each item's links sit next to the item name in the first column and open each figure (`png`, `pdf`) or table (`md`, a rendered table; the data file in its own format where no markdown rendering is committed) directly from the repository; the number registry and the IxIDN-Matched quality-control row link their `tsv` files. They are ordinary relative links, so they work on GitHub and in a local checkout; Ctrl/Cmd-click or middle-click opens one in a new tab (plain Markdown cannot force a new tab). Figures S4, S6 and S26 are committed as `png` only. A row with no link is not a single figure or table and has no committed `tsv` (the SI numbers are in `*_numbers.json`). Script and Output entries show the file name and link to the file; the full relative path is the link target.

### Main text

<table>
<thead>
<tr><th width="14%">Item</th><th width="52%">Content</th><th width="16%">Script</th><th width="18%">Output</th></tr>
</thead>
<tbody>
<tr><td>Figure 1 <a href="figures/main/fig1_knowledge_transfer.png">[png]</a> <a href="figures/main/fig1_knowledge_transfer.pdf">[pdf]</a></td><td>Knowledge transfer between rare and common diseases (panels a, b from Open Targets 26.09; c, d schematic)</td><td><a href="figures/main/fig1_knowledge_transfer.py"><code>fig1_<wbr>knowledge_<wbr>transfer.py</code></a>, drawing code <a href="src/disatlas/fig1_open_targets.py"><code>fig1_<wbr>open_<wbr>targets.py</code></a>; reads the committed Open Targets 26.09 tables <a href="data/derived/fig1_open_targets/"><code>fig1_<wbr>open_<wbr>targets/</code></a> (rebuilding them from the release is the module's <code>tables</code> subcommand and is not part of <code>run_pipeline.sh</code>)</td><td><a href="figures/main/fig1_knowledge_transfer.png"><code>fig1_<wbr>knowledge_<wbr>transfer.png</code></a> <a href="figures/main/fig1_knowledge_transfer.pdf"><code>.pdf</code></a>; the same script also writes Figure S1</td></tr>
<tr><td>Figure 2 <a href="figures/main/disatlas_pipeline.png">[png]</a> <a href="figures/main/disatlas_pipeline.pdf">[pdf]</a></td><td>DisAtlas pipeline and its two evaluations</td><td><a href="figures/main/fig2_disatlas_pipeline.py"><code>fig2_<wbr>disatlas_<wbr>pipeline.py</code></a></td><td><a href="figures/main/disatlas_pipeline.png"><code>disatlas_<wbr>pipeline.png</code></a> <a href="figures/main/disatlas_pipeline.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure 3 <a href="figures/main/ixidn_benchmark.png">[png]</a> <a href="figures/main/ixidn_benchmark.pdf">[pdf]</a></td><td>The IxIDN-Matched benchmark</td><td><a href="figures/main/fig3_ixidn_benchmark.py"><code>fig3_<wbr>ixidn_<wbr>benchmark.py</code></a></td><td><a href="figures/main/ixidn_benchmark.png"><code>ixidn_<wbr>benchmark.png</code></a> <a href="figures/main/ixidn_benchmark.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure 4 <a href="figures/main/fig4_evidence_trace.png">[png]</a> <a href="figures/main/fig4_evidence_trace.pdf">[pdf]</a></td><td>Tracing a disease-disease similarity to its evidence</td><td><a href="figures/main/fig4_evidence_trace.py"><code>fig4_<wbr>evidence_<wbr>trace.py</code></a>; prose numbers <a href="figures/main/fig4_numbers.py"><code>fig4_<wbr>numbers.py</code></a></td><td><a href="figures/main/fig4_evidence_trace.png"><code>fig4_<wbr>evidence_<wbr>trace.png</code></a> <a href="figures/main/fig4_evidence_trace.pdf"><code>.pdf</code></a>, <a href="figures/main/fig4_numbers.json"><code>fig4_<wbr>numbers.json</code></a></td></tr>
<tr><td>Figure 5 <a href="figures/main/disease_clustering_map.png">[png]</a> <a href="figures/main/disease_clustering_map.pdf">[pdf]</a></td><td>Disease map, clusters, and DisAtlas vs TxGNN by disease rarity (panel c)</td><td><a href="figures/main/fig5_disease_clustering_map.py"><code>fig5_<wbr>disease_<wbr>clustering_<wbr>map.py</code></a></td><td><a href="figures/main/disease_clustering_map.png"><code>disease_<wbr>clustering_<wbr>map.png</code></a> <a href="figures/main/disease_clustering_map.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure 6 <a href="figures/main/fig5_method_comparison.png">[png]</a> <a href="figures/main/fig5_method_comparison.pdf">[pdf]</a></td><td>Cross-system comparison on IxIDN-Matched</td><td><a href="figures/main/fig6_method_comparison.py"><code>fig6_<wbr>method_<wbr>comparison.py</code></a></td><td><a href="figures/main/fig5_method_comparison.png"><code>fig5_<wbr>method_<wbr>comparison.png</code></a> <a href="figures/main/fig5_method_comparison.pdf"><code>.pdf</code></a> (same file name as the manuscript figure)</td></tr>
<tr><td>Figure 7 <a href="figures/main/fig7_prospective.png">[png]</a> <a href="figures/main/fig7_prospective.pdf">[pdf]</a></td><td>Event-time test: rank of the case-study disease and conditional-logit odds ratios by year and method</td><td><a href="figures/main/fig7_prospective.py"><code>fig7_<wbr>prospective.py</code></a></td><td><a href="figures/main/fig7_prospective.png"><code>fig7_<wbr>prospective.png</code></a> <a href="figures/main/fig7_prospective.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Table 1 <a href="tables/main/table1_protocols.md">[md]</a></td><td>Protocols of the cross-system comparison</td><td><a href="tables/main/table1_protocols.py"><code>table1_<wbr>protocols.py</code></a></td><td><a href="tables/main/table1_protocols.md"><code>table1_<wbr>protocols.md</code></a> <a href="tables/main/table1_protocols.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Main-text numbers <a href="tests/main_text_numbers.tsv">[tsv]</a></td><td>All 568 registered number tokens of the main text, each with its kind and check</td><td><a href="tests/main_quantities.py"><code>main_<wbr>quantities.py</code></a>, <a href="tests/test_main_parity.py"><code>test_<wbr>main_<wbr>parity.py</code></a></td><td>registry <a href="tests/main_text_numbers.tsv"><code>main_<wbr>text_<wbr>numbers.tsv</code></a></td></tr>
</tbody>
</table>

### Supplementary figures

<table>
<thead>
<tr><th width="14%">Item</th><th width="52%">Content</th><th width="16%">Script</th><th width="18%">Output</th></tr>
</thead>
<tbody>
<tr><td>Figure S1 <a href="figures/SI/rare_disease_stats.png">[png]</a> <a href="figures/SI/rare_disease_stats.pdf">[pdf]</a></td><td>Prevalence of rare diseases and sources of their genetic evidence</td><td><a href="figures/main/fig1_knowledge_transfer.py"><code>fig1_<wbr>knowledge_<wbr>transfer.py</code></a> (panel a is the literature prevalence panel, used unchanged)</td><td><a href="figures/SI/rare_disease_stats.png"><code>rare_<wbr>disease_<wbr>stats.png</code></a> <a href="figures/SI/rare_disease_stats.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S2 <a href="figures/SI/SI_approved_repurposing_schematic.png">[png]</a> <a href="figures/SI/SI_approved_repurposing_schematic.pdf">[pdf]</a></td><td>Approved-repurposing benchmark (schematic)</td><td><a href="figures/SI/si_fda_schematic.py"><code>si_<wbr>fda_<wbr>schematic.py</code></a></td><td><a href="figures/SI/SI_approved_repurposing_schematic.png"><code>SI_<wbr>approved_<wbr>repurposing_<wbr>schematic.png</code></a> <a href="figures/SI/SI_approved_repurposing_schematic.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S3 <a href="figures/SI/SI_rare_common_transfer.png">[png]</a> <a href="figures/SI/SI_rare_common_transfer.pdf">[pdf]</a></td><td>DisAtlas performance by disease rarity and by anchor disease</td><td><a href="figures/SI/si_rare_common_transfer.py"><code>si_<wbr>rare_<wbr>common_<wbr>transfer.py</code></a></td><td><a href="figures/SI/SI_rare_common_transfer.png"><code>SI_<wbr>rare_<wbr>common_<wbr>transfer.png</code></a> <a href="figures/SI/SI_rare_common_transfer.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S4 <a href="figures/SI/LCA.png">[png]</a></td><td>Distance to the lowest common ancestor (illustration)</td><td>none (static image)</td><td><a href="figures/SI/LCA.png"><code>LCA.png</code></a></td></tr>
<tr><td>Figure S5 <a href="figures/SI/SI_graph_channels_schematic.png">[png]</a> <a href="figures/SI/SI_graph_channels_schematic.pdf">[pdf]</a></td><td>How the DisAtlas training graph is built and scored</td><td><a href="figures/SI/si_graph_channels_schematic.py"><code>si_<wbr>graph_<wbr>channels_<wbr>schematic.py</code></a></td><td><a href="figures/SI/SI_graph_channels_schematic.png"><code>SI_<wbr>graph_<wbr>channels_<wbr>schematic.png</code></a> <a href="figures/SI/SI_graph_channels_schematic.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S6 <a href="figures/SI/SI_mechdist_embedding_pair_score.png">[png]</a></td><td>From the MechDist network to embeddings and pair scores; separation on IxIDN-Matched</td><td><a href="figures/SI/si_mechdist_embedding_pair_score.py"><code>si_<wbr>mechdist_<wbr>embedding_<wbr>pair_<wbr>score.py</code></a></td><td><a href="figures/SI/SI_mechdist_embedding_pair_score.png"><code>SI_<wbr>mechdist_<wbr>embedding_<wbr>pair_<wbr>score.png</code></a>, <a href="figures/SI/SI_mechdist_embedding_pair_score_values.json"><code>_<wbr>values.json</code></a></td></tr>
<tr><td>Figure S7 <a href="figures/SI/si_ixidn_matched_funnel.png">[png]</a> <a href="figures/SI/si_ixidn_matched_funnel.pdf">[pdf]</a></td><td>Construction of IxIDN-Matched and its shortcut audit</td><td><a href="figures/SI/si_ixidn_matched_funnel.py"><code>si_<wbr>ixidn_<wbr>matched_<wbr>funnel.py</code></a></td><td><a href="figures/SI/si_ixidn_matched_funnel.png"><code>si_<wbr>ixidn_<wbr>matched_<wbr>funnel.png</code></a> <a href="figures/SI/si_ixidn_matched_funnel.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S8 <a href="figures/SI/si_pipeline_detailed.png">[png]</a> <a href="figures/SI/si_pipeline_detailed.pdf">[pdf]</a></td><td>Detailed DisAtlas pipeline</td><td><a href="figures/SI/si_pipeline_detailed.py"><code>si_<wbr>pipeline_<wbr>detailed.py</code></a></td><td><a href="figures/SI/si_pipeline_detailed.png"><code>si_<wbr>pipeline_<wbr>detailed.png</code></a> <a href="figures/SI/si_pipeline_detailed.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S9 <a href="figures/SI/si_mechanism_sensitivity.png">[png]</a> <a href="figures/SI/si_mechanism_sensitivity.pdf">[pdf]</a></td><td>Broad and specific evidence channels of the MechDist network</td><td><a href="figures/SI/si_mechanism_sensitivity.py"><code>si_<wbr>mechanism_<wbr>sensitivity.py</code></a></td><td><a href="figures/SI/si_mechanism_sensitivity.png"><code>si_<wbr>mechanism_<wbr>sensitivity.png</code></a> <a href="figures/SI/si_mechanism_sensitivity.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S10 <a href="figures/SI/si_mechdist_degree_by_channel.png">[png]</a> <a href="figures/SI/si_mechdist_degree_by_channel.pdf">[pdf]</a></td><td>Degree in each graph of the MechDist network</td><td><a href="figures/SI/si_mechdist_degree_by_channel.py"><code>si_<wbr>mechdist_<wbr>degree_<wbr>by_<wbr>channel.py</code></a></td><td><a href="figures/SI/si_mechdist_degree_by_channel.png"><code>si_<wbr>mechdist_<wbr>degree_<wbr>by_<wbr>channel.png</code></a> <a href="figures/SI/si_mechdist_degree_by_channel.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S11 <a href="figures/SI/si_threshold_sensitivity.png">[png]</a> <a href="figures/SI/si_threshold_sensitivity.pdf">[pdf]</a></td><td>Threshold sensitivity</td><td><a href="figures/SI/si_threshold_sensitivity.py"><code>si_<wbr>threshold_<wbr>sensitivity.py</code></a></td><td><a href="figures/SI/si_threshold_sensitivity.png"><code>si_<wbr>threshold_<wbr>sensitivity.png</code></a> <a href="figures/SI/si_threshold_sensitivity.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S12 <a href="figures/SI/si_ablation_31arm.png">[png]</a> <a href="figures/SI/si_ablation_31arm.pdf">[pdf]</a></td><td>Component ablations (31 arms)</td><td><a href="figures/SI/si_ablation_31arm.py"><code>si_<wbr>ablation_<wbr>31arm.py</code></a></td><td><a href="figures/SI/si_ablation_31arm.png"><code>si_<wbr>ablation_<wbr>31arm.png</code></a> <a href="figures/SI/si_ablation_31arm.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S13 <a href="figures/SI/si_uncertainty.png">[png]</a> <a href="figures/SI/si_uncertainty.pdf">[pdf]</a></td><td>Seed and sampling uncertainty</td><td><a href="figures/SI/si_uncertainty.py"><code>si_<wbr>uncertainty.py</code></a></td><td><a href="figures/SI/si_uncertainty.png"><code>si_<wbr>uncertainty.png</code></a> <a href="figures/SI/si_uncertainty.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S14 <a href="figures/SI/si_placement_confidence.png">[png]</a> <a href="figures/SI/si_placement_confidence.pdf">[pdf]</a></td><td>Per-disease placement confidence: exclusion sensitivity</td><td><a href="figures/SI/si_placement_confidence.py"><code>si_<wbr>placement_<wbr>confidence.py</code></a></td><td><a href="figures/SI/si_placement_confidence.png"><code>si_<wbr>placement_<wbr>confidence.png</code></a> <a href="figures/SI/si_placement_confidence.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S15 <a href="figures/SI/si_node_disjoint.png">[png]</a> <a href="figures/SI/si_node_disjoint.pdf">[pdf]</a></td><td>Node-disjoint evaluation</td><td><a href="figures/SI/si_node_disjoint.py"><code>si_<wbr>node_<wbr>disjoint.py</code></a></td><td><a href="figures/SI/si_node_disjoint.png"><code>si_<wbr>node_<wbr>disjoint.png</code></a> <a href="figures/SI/si_node_disjoint.pdf"><code>.pdf</code></a>, <a href="figures/SI/si_node_disjoint_values.tsv"><code>si_<wbr>node_<wbr>disjoint_<wbr>values.tsv</code></a> <a href="figures/SI/si_node_disjoint_values.json"><code>.json</code></a></td></tr>
<tr><td>Figure S16 <a href="figures/SI/extended_cluster_validation.png">[png]</a> <a href="figures/SI/extended_cluster_validation.pdf">[pdf]</a></td><td>Extended cluster validation</td><td><a href="figures/SI/si_cluster_validation.py"><code>si_<wbr>cluster_<wbr>validation.py</code></a></td><td><a href="figures/SI/extended_cluster_validation.png"><code>extended_<wbr>cluster_<wbr>validation.png</code></a> <a href="figures/SI/extended_cluster_validation.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S17 <a href="figures/SI/SI_cluster_map_vs_model_distance.png">[png]</a> <a href="figures/SI/SI_cluster_map_vs_model_distance.pdf">[pdf]</a></td><td>Display layout fidelity</td><td><a href="figures/SI/si_cluster_map_vs_model.py"><code>si_<wbr>cluster_<wbr>map_<wbr>vs_<wbr>model.py</code></a></td><td><a href="figures/SI/SI_cluster_map_vs_model_distance.png"><code>SI_<wbr>cluster_<wbr>map_<wbr>vs_<wbr>model_<wbr>distance.png</code></a> <a href="figures/SI/SI_cluster_map_vs_model_distance.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S18 <a href="figures/SI/SI_rare_common_neighbour_coverage.png">[png]</a> <a href="figures/SI/SI_rare_common_neighbour_coverage.pdf">[pdf]</a></td><td>Rare–common neighbor coverage</td><td><a href="figures/SI/si_rare_common_neighbour_coverage.py"><code>si_<wbr>rare_<wbr>common_<wbr>neighbour_<wbr>coverage.py</code></a></td><td><a href="figures/SI/SI_rare_common_neighbour_coverage.png"><code>SI_<wbr>rare_<wbr>common_<wbr>neighbour_<wbr>coverage.png</code></a> <a href="figures/SI/SI_rare_common_neighbour_coverage.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S19 <a href="figures/SI/SI_event_time_auc_methods.png">[png]</a> <a href="figures/SI/SI_event_time_auc_methods.pdf">[pdf]</a></td><td>Event-time holdout: per-drug probability by method and decision year (matched analysis)</td><td><a href="figures/SI/si_event_time_auc_methods.py"><code>si_<wbr>event_<wbr>time_<wbr>auc_<wbr>methods.py</code></a></td><td><a href="figures/SI/SI_event_time_auc_methods.png"><code>SI_<wbr>event_<wbr>time_<wbr>auc_<wbr>methods.png</code></a> <a href="figures/SI/SI_event_time_auc_methods.pdf"><code>.pdf</code></a>, <a href="figures/SI/SI_event_time_auc_methods_values.json"><code>SI_<wbr>event_<wbr>time_<wbr>auc_<wbr>methods_<wbr>values.json</code></a></td></tr>
<tr><td>Figure S20 <a href="figures/SI/SI_event_time_hits10.png">[png]</a> <a href="figures/SI/SI_event_time_hits10.pdf">[pdf]</a></td><td>Event-time holdout: secondary endpoint with and without matching</td><td><a href="figures/main/fig7_prospective.py"><code>fig7_<wbr>prospective.py</code></a> <code>--si hits10</code></td><td><a href="figures/SI/SI_event_time_hits10.png"><code>SI_<wbr>event_<wbr>time_<wbr>hits10.png</code></a> <a href="figures/SI/SI_event_time_hits10.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S21 <a href="figures/SI/SI_event_time_weighted.png">[png]</a> <a href="figures/SI/SI_event_time_weighted.pdf">[pdf]</a></td><td>Event-time holdout: weighted for disease popularity (post hoc)</td><td><a href="figures/main/fig7_prospective.py"><code>fig7_<wbr>prospective.py</code></a> <code>--si weighted</code></td><td><a href="figures/SI/SI_event_time_weighted.png"><code>SI_<wbr>event_<wbr>time_<wbr>weighted.png</code></a> <a href="figures/SI/SI_event_time_weighted.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S22 <a href="figures/SI/si_popularity_strata.png">[png]</a> <a href="figures/SI/si_popularity_strata.pdf">[pdf]</a></td><td>Later-trialed diseases within popularity strata</td><td><a href="figures/SI/si_popularity_strata.py"><code>si_<wbr>popularity_<wbr>strata.py</code></a></td><td><a href="figures/SI/si_popularity_strata.png"><code>si_<wbr>popularity_<wbr>strata.png</code></a> <a href="figures/SI/si_popularity_strata.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S23 <a href="figures/SI/SI_fda_channels.png">[png]</a> <a href="figures/SI/SI_fda_channels.pdf">[pdf]</a></td><td>Broad and specific evidence channels on FDA-approved repurposing pairs</td><td><a href="figures/SI/si_fda_channels.py"><code>si_<wbr>fda_<wbr>channels.py</code></a></td><td><a href="figures/SI/SI_fda_channels.png"><code>SI_<wbr>fda_<wbr>channels.png</code></a> <a href="figures/SI/SI_fda_channels.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S24 <a href="figures/SI/SI_complement_milestones.png">[png]</a> <a href="figures/SI/SI_complement_milestones.pdf">[pdf]</a></td><td>Clinical milestones of complement inhibitors and DisAtlas similarity of the eculizumab indications</td><td><a href="figures/SI/si_complement_milestones.py"><code>si_<wbr>complement_<wbr>milestones.py</code></a></td><td><a href="figures/SI/SI_complement_milestones.png"><code>SI_<wbr>complement_<wbr>milestones.png</code></a> <a href="figures/SI/SI_complement_milestones.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S25 <a href="figures/SI/SI_mechanism_coverage_slices.png">[png]</a> <a href="figures/SI/SI_mechanism_coverage_slices.pdf">[pdf]</a></td><td>Mechanism coverage of approved expansions in dated graph slices</td><td><a href="figures/SI/si_mechanism_coverage_slices.py"><code>si_<wbr>mechanism_<wbr>coverage_<wbr>slices.py</code></a></td><td><a href="figures/SI/SI_mechanism_coverage_slices.png"><code>SI_<wbr>mechanism_<wbr>coverage_<wbr>slices.png</code></a> <a href="figures/SI/SI_mechanism_coverage_slices.pdf"><code>.pdf</code></a></td></tr>
<tr><td>Figure S26 <a href="figures/SI/SI_explainability.png">[png]</a></td><td>Evidence behind three rare–common case pairs (cosine scores; TxGNN panels static)</td><td><a href="figures/SI/si_explainability.py"><code>si_<wbr>explainability.py</code></a></td><td><a href="figures/SI/SI_explainability.png"><code>SI_<wbr>explainability.png</code></a>, <a href="figures/SI/SI_explainability_values.json"><code>SI_<wbr>explainability_<wbr>values.json</code></a></td></tr>
<tr><td>Methods, Figure 3b, Figure S7 <a href="data/raw/benchmarks/ixidn_matched/llm_curation_verdicts.tsv" title="llm_curation_verdicts.tsv">[tsv]</a> <a href="data/derived/fig3_ixidn/qc_verdict_summary.tsv" title="qc_verdict_summary.tsv">[tsv]</a></td><td>IxIDN-Matched quality control: language-model verdicts on the 4,302 candidate positives (1,285 strong + 897 plausible kept; repeat check 84.7%) and the seeded cell match</td><td><a href="src/disatlas/benchmarks/ixidn_matched.py"><code>ixidn_<wbr>matched.py</code></a> (stored verdicts <a href="data/raw/benchmarks/ixidn_matched/llm_curation_verdicts.tsv"><code>llm_<wbr>curation_<wbr>verdicts.tsv</code></a>), rerun <a href="src/disatlas/benchmarks/matched_qc.py"><code>matched_<wbr>qc.py</code></a>, <a href="tests/test_ixidn_matched_qc.py"><code>test_<wbr>ixidn_<wbr>matched_<wbr>qc.py</code></a></td><td><a href="data/derived/fig3_ixidn/qc_verdict_summary.tsv"><code>qc_<wbr>verdict_<wbr>summary.tsv</code></a> (kept/removed by reason)</td></tr>
</tbody>
</table>

### Supplementary tables

<table>
<thead>
<tr><th width="14%">Item</th><th width="52%">Content</th><th width="16%">Script</th><th width="18%">Output</th></tr>
</thead>
<tbody>
<tr><td>Table S1 <a href="tables/SI/table_input_dis_gene_sources.md">[md]</a></td><td>Disease-gene association sources</td><td>none (descriptive)</td><td><a href="tables/SI/table_input_dis_gene_sources.md"><code>table_<wbr>input_<wbr>dis_<wbr>gene_<wbr>sources.md</code></a> <a href="tables/SI/table_input_dis_gene_sources.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S2 <a href="tables/SI/table_input_dis_symptom_sources.md">[md]</a></td><td>Disease-symptom association sources</td><td>none (descriptive)</td><td><a href="tables/SI/table_input_dis_symptom_sources.md"><code>table_<wbr>input_<wbr>dis_<wbr>symptom_<wbr>sources.md</code></a> <a href="tables/SI/table_input_dis_symptom_sources.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S3 <a href="tables/SI/table_fig3a_trials.md">[md]</a></td><td>Trials behind Figure 3a</td><td><a href="tables/SI/table_fig3a_trials.py"><code>table_<wbr>fig3a_<wbr>trials.py</code></a></td><td><a href="tables/SI/table_fig3a_trials.md"><code>table_<wbr>fig3a_<wbr>trials.md</code></a> <a href="tables/SI/table_fig3a_trials.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S4 <a href="tables/SI/table_ontology_baseline.md">[md]</a></td><td>Ontology-only baselines and concordance checks</td><td><a href="tables/SI/table_ontology_baseline.py"><code>table_<wbr>ontology_<wbr>baseline.py</code></a></td><td><a href="tables/SI/table_ontology_baseline.md"><code>table_<wbr>ontology_<wbr>baseline.md</code></a></td></tr>
<tr><td>Table S5 <a href="tables/SI/table_si_baselines.md">[md]</a></td><td>Contribution of the network and of its representation</td><td><a href="tables/SI/table_si_baselines.py"><code>table_<wbr>si_<wbr>baselines.py</code></a></td><td><a href="tables/SI/table_si_baselines.md"><code>table_<wbr>si_<wbr>baselines.md</code></a> <a href="tables/SI/table_si_baselines.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S6 <a href="tables/SI/table_llm_subgroups.md">[md]</a></td><td>Node-disjoint DisAtlas vs language-model systems by subgroup</td><td><a href="tables/SI/table_llm_subgroups.py"><code>table_<wbr>llm_<wbr>subgroups.py</code></a></td><td><a href="tables/SI/table_llm_subgroups.md"><code>table_<wbr>llm_<wbr>subgroups.md</code></a> <a href="tables/SI/table_llm_subgroups.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S7 <a href="tables/SI/table_heuristic_vs_ai.md">[md]</a></td><td>Score computed on the evidence network vs each AI system</td><td><a href="tables/SI/table_heuristic_vs_ai.py"><code>table_<wbr>heuristic_<wbr>vs_<wbr>ai.py</code></a></td><td><a href="tables/SI/table_heuristic_vs_ai.md"><code>table_<wbr>heuristic_<wbr>vs_<wbr>ai.md</code></a> <a href="tables/SI/table_heuristic_vs_ai.tsv"><code>.tsv</code></a>, <a href="tables/SI/table_heuristic_vs_ai_numbers.json"><code>_<wbr>numbers.json</code></a></td></tr>
<tr><td>Table S8 <a href="tables/SI/table_degree_controls.md">[md]</a></td><td>Degree-matched and within-disease AUROC</td><td><a href="tables/SI/table_degree_controls.py"><code>table_<wbr>degree_<wbr>controls.py</code></a></td><td><a href="tables/SI/table_degree_controls.md"><code>table_<wbr>degree_<wbr>controls.md</code></a> <a href="tables/SI/table_degree_controls.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S9 <a href="tables/SI/table_cluster_controls.md">[md]</a></td><td>Positive and negative control pairs</td><td><a href="tables/SI/table_cluster_controls.py"><code>table_<wbr>cluster_<wbr>controls.py</code></a></td><td><a href="tables/SI/table_cluster_controls.md"><code>table_<wbr>cluster_<wbr>controls.md</code></a> <a href="tables/SI/table_cluster_controls.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S10 <a href="tables/SI/table_cluster_themes.md">[md]</a></td><td>Theme of each of the 55 disease clusters</td><td><a href="tables/SI/table_cluster_themes.py"><code>table_<wbr>cluster_<wbr>themes.py</code></a></td><td><a href="tables/SI/table_cluster_themes.md"><code>table_<wbr>cluster_<wbr>themes.md</code></a> <a href="tables/SI/table_cluster_themes.tsv"><code>.tsv</code></a>, <a href="tables/SI/table_cluster_themes_numbers.json"><code>_<wbr>numbers.json</code></a></td></tr>
<tr><td>Table S11 <a href="tables/SI/table_rare_common_strata.md">[md]</a></td><td>AUROC by disease-rarity pairing</td><td><a href="tables/SI/table_rare_common_strata.py"><code>table_<wbr>rare_<wbr>common_<wbr>strata.py</code></a></td><td><a href="tables/SI/table_rare_common_strata.md"><code>table_<wbr>rare_<wbr>common_<wbr>strata.md</code></a> <a href="tables/SI/table_rare_common_strata.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S12 <a href="tables/SI/table_ra_myositis.md">[md]</a></td><td>RA against the inflammatory myopathies and overlap disorders</td><td><a href="tables/SI/table_ra_myositis.py"><code>table_<wbr>ra_<wbr>myositis.py</code></a> (values from <a href="figures/SI/si_cluster_rarecommon_numbers.py"><code>si_<wbr>cluster_<wbr>rarecommon_<wbr>numbers.py</code></a>)</td><td><a href="tables/SI/table_ra_myositis.md"><code>table_<wbr>ra_<wbr>myositis.md</code></a> <a href="tables/SI/table_ra_myositis.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Tables S13, S14 <a href="figures/SI/SI_event_time_tables.md">[md]</a></td><td>Event-time holdout: pooled advantages; subsets of new trials</td><td><a href="figures/SI/si_event_time_tables.py"><code>si_<wbr>event_<wbr>time_<wbr>tables.py</code></a></td><td><a href="figures/SI/SI_event_time_tables.md"><code>SI_<wbr>event_<wbr>time_<wbr>tables.md</code></a></td></tr>
<tr><td>Tables S15, S16 <a href="figures/SI/SI_event_time_or_tables.md">[md]</a></td><td>Event-time: eculizumab rank of MS and NMO by year; pooled per-drug probability scale</td><td><a href="figures/SI/si_event_time_or_tables.py"><code>si_<wbr>event_<wbr>time_<wbr>or_<wbr>tables.py</code></a></td><td><a href="figures/SI/SI_event_time_or_tables.md"><code>SI_<wbr>event_<wbr>time_<wbr>or_<wbr>tables.md</code></a>, <a href="figures/SI/SI_event_time_or_tables_numbers.json"><code>_<wbr>numbers.json</code></a></td></tr>
<tr><td>Table S17 <a href="figures/SI/SI_go_release_table.md">[md]</a></td><td>Advantage with current and slice-dated GO annotations</td><td><a href="figures/SI/si_go_release_table.py"><code>si_<wbr>go_<wbr>release_<wbr>table.py</code></a></td><td><a href="figures/SI/SI_go_release_table.md"><code>SI_<wbr>go_<wbr>release_<wbr>table.md</code></a> <a href="figures/SI/SI_go_release_table.tsv"><code>.tsv</code></a></td></tr>
<tr><td>Table S18 <a href="tables/SI/table_SI_trace_order.md">[md]</a></td><td>Figure 4 intermediates changed by specificity weighting</td><td><a href="tables/SI/table_SI_trace_order.py"><code>table_<wbr>SI_<wbr>trace_<wbr>order.py</code></a></td><td><a href="tables/SI/table_SI_trace_order.md"><code>table_<wbr>SI_<wbr>trace_<wbr>order.md</code></a> <a href="tables/SI/table_SI_trace_order.tsv"><code>.tsv</code></a> (displayed-item values: <a href="tables/SI/table_SI_trace_order_values.md"><code>table_<wbr>SI_<wbr>trace_<wbr>order_<wbr>values.md</code></a> <a href="tables/SI/table_SI_trace_order_values.tsv"><code>.tsv</code></a>)</td></tr>
<tr><td>SI numbers</td><td>Numbers quoted in the SI text</td><td><code>figures/SI/si_<wbr>*_<wbr>numbers.py</code></td><td><code>figures/SI/SI_<wbr>*_<wbr>numbers.json</code>; asserted by <a href="tests/test_sibench_parity.py"><code>test_<wbr>sibench_<wbr>parity.py</code></a>, <a href="tests/test_sibio_parity.py"><code>test_<wbr>sibio_<wbr>parity.py</code></a> (specs in <a href="tests/parity/"><code>parity/</code></a>)</td></tr>
</tbody>
</table>

## Repository structure

```
run_pipeline.sh          one command for figures, tables, tests and notebooks (see Run modes)
requirements.txt         pinned Python packages (warm pipeline; no torch)
requirements-retrain.txt torch, only for the optional cold retrain
data/
  raw/                   inputs of the graph build: public disease-gene sources, Open Targets 26.09, HPO,
                         the 8,343 released knowledge-graph disease-gene edges, broad-channel gene tables and the
                         released specific-channel edges (channels/specific/); see data/raw/README.md
  raw/benchmarks/        IxIDN+/IxIDN- pairs, IxIDN-Matched pairs and curation verdicts, FDA controls, comparator outputs
  derived/               small precomputed result tables read by the figure, table and test scripts
  zenodo/                (download, not in git) broad graph, term-distance slices, frozen model and its outputs,
                         score matrices, large tables
  README.md              column dictionaries, licence and blinding
  ZENODO_SHA256SUMS.txt  checksums of every data/zenodo file (manifest with sizes and licences: release/ZENODO_MANIFEST.tsv)
src/disatlas/            library: graph build, model loading and scoring, retraining, blinding, paths, figure style
figures/main, figures/SI scripts for the main and Supplementary figures, and their outputs
tables/main, tables/SI   scripts for Table 1 and the Supplementary tables, and their outputs
notebooks/               four exploratory notebooks (data and graph, IxIDN benchmarks, SI robustness, SI biology)
tests/                   parity tests and the main-text number registry
scripts/                 P_full parity replay, notebook runner; scripts/maintainer/ shows how the
                         shipped data were exported (not needed to reproduce)
docs/                    model outputs, optional P_full retrain, errata
release/                 data-minimisation and knowledge-graph exposure ledgers
```

### Training graph

The training graph (P_full: 532,958 edges, 13,042 diseases) is the union of two disease-disease channels:

- **broad** (369,405 edges, shipped as `data/zenodo/graph/graph_edges.parquet`): its disease-gene table of record
  has 24,450 tier-1 rows, 604 of them licence-restricted pairs that are not redistributed (see
  [Data and licensing](#data-and-licensing)). `python -m disatlas.graph.raw_build` rebuilds the broad channel from
  `data/raw`, whose disease-gene table is 21,066 rows from public sources (GenCC, Orphanet, Open Targets 26.09
  direct) plus the 8,343 released knowledge-graph disease-gene edges, and nothing else. Of its 23,840 tier-1 rows,
  16,975 are public and 6,865 come from the released edges; the rebuilt broad graph has 364,288 edges. The GO and HPO term-distance slices it reads ship in
  `data/zenodo/ontology_term_term_distance/`.
- **specific** (189,903 edges): shipped as disease-disease edges with their weight only
  (`data/raw/channels/specific/specific_edges.parquet`). Rebuilding it from its inputs would need about 1.06 million
  rows of licence-restricted knowledge-graph evidence (literature, pathology, MetaBase, CTD and mechanism nodes),
  which is why only its edges are released; a rebuild from public sources only is FAR from the released model
  (test AUROC -0.066).

`python -m disatlas.graph.pfull_build` combines them (outer union, FDA-pair strip, walk weight `w_arm`,
gene-anchored step) and reproduces the recorded P_full graph exactly (edge set, row order and `w_arm`,
`tests/test_pfull_graph.py`). The released disease-gene edges alone are not sufficient: the training graph also
needs the specific-channel edges.

## Quick start

Requirements: Linux or macOS, Python 3.11, about 16 GB of RAM, and about 1.2 GB of disk for the unpacked Zenodo bundle.

```bash
git clone <repository URL> DisAtlas && cd DisAtlas
python3.11 -m venv .venv
.venv/bin/pip install -r requirements.txt
```

`requirements.txt` pins the package versions used to render the figures; with them a fresh run regenerates every committed
figure byte for byte (other matplotlib versions can shift anti-aliasing by a few grey levels without changing any number).

Download the Zenodo record (DOI: [10.5281/zenodo.24263240](https://doi.org/10.5281/zenodo.24263240)) and unpack it so that
`data/zenodo/` contains `graph/`, `ontology_term_term_distance/`, `model/`, `model_outputs/`, `figures_main/`, `si_bio/`, `event_time/`. Then check it:

```bash
cd data/zenodo && sha256sum -c ../ZENODO_SHA256SUMS.txt && cd ../..   # macOS: shasum -a 256 -c
./run_pipeline.sh
```

`run_pipeline.sh` uses `.venv/bin/python` (override with `PYTHON=/path/to/python`) and sets `PYTHONPATH=src`.
To run a single script by hand: `export PYTHONPATH=$PWD/src MPLBACKEND=Agg` and then, for example,
`.venv/bin/python figures/main/fig3_ixidn_benchmark.py`.

## Run modes

| Mode | Command | What it does | Agreement with the manuscript |
|---|---|---|---|
| WARM (default) | `./run_pipeline.sh` | Uses the processed data and frozen model outputs from `data/zenodo`; rebuilds the P_full training graph and IxIDN-Matched as checks | exact |
| COLD | `./run_pipeline.sh --cold` | Rebuilds the P_full training graph into `build/` (with `DISATLAS_RAWGRAPH=1` also the broad graph from `data/raw`), retrains the 5-seed ensemble and re-scores the benchmarks, checks the graph against its recorded signature and the retrain with the pre-registered parity test, then makes all figures and tables and runs the tests | graph exact; retrain judged by the parity test |
| QC rerun (opt-in) | `./run_pipeline.sh --cold-qc` (either mode) | Stage `qc` always rebuilds IxIDN-Matched from the stored QC verdicts and checks it equals the shipped file; `--cold-qc` also reruns the QC with your own key (`ANTHROPIC_API_KEY`, `QC_MODEL`) into `build/qc/` using a prompt reconstructed from the Methods; drug names and trial phases from the licensed trial source are optional inputs (`--drugs`, `--phases`) | stored verdicts exact; rerun at a tolerance (`tests/test_ixidn_matched_qc.py`) |

Figures, tables, notebooks and the manuscript numbers always come from the frozen record (they need
`data/zenodo`), in both modes, because a retrained ensemble reproduces the published scores closely but not bit for bit.

Stages, in order:

| Mode | Stages |
|---|---|
| WARM | `verify graph qc figures si tests notebooks` |
| COLD | `verify graph qc retrain score figures si tests notebooks` |

- `verify`: checks `data/zenodo` against `data/ZENODO_SHA256SUMS.txt` (in COLD mode it is skipped with a message if `data/zenodo` is absent; comparisons with the frozen record are then skipped).
- `graph`: `python -m disatlas.graph.pfull_build --out build/graph` (broad graph from `data/zenodo/graph`; in COLD
  mode with `DISATLAS_RAWGRAPH=1`, first `python -m disatlas.graph.raw_build --out build/rawgraph`, about 1 h, and
  the P_full graph is built from that). Both builds check their output against the recorded signatures.
- `qc`: `python -m disatlas.benchmarks.ixidn_matched --out build/ixidn_matched --check` (stored verdicts); with `--cold-qc` also `python -m disatlas.benchmarks.matched_qc --out build/qc --model $QC_MODEL`.
- `retrain` (COLD): `python -m disatlas.retrain --graph build/graph/pfull_graph_edges.parquet`, one process per seed (4707-4711) in parallel; weighted node2vec on `w_arm`, dim 32, walk length 20, window 2, p = 1, q = 0.5; pairs are scored by the 5-seed mean embedding cosine (no decoder). Needs `requirements-retrain.txt` (torch) on top of `requirements.txt`.
- `score` (COLD): `python -m disatlas.rescore --out build/scores --all-pairs`.
- `figures`, `si`: all scripts in `figures/main`, `tables/main`, then `figures/SI`, `tables/SI`.
- `tests`: `pytest tests` (in COLD mode first `tests/test_cold_exact.py` and `tests/test_cold_tolerance.py` on `build/`).
- `notebooks`: executes the four notebooks.

Start at a stage with `--from STAGE` and run one stage with `--only STAGE`, for example
`./run_pipeline.sh --only tests` or `./run_pipeline.sh --cold --from retrain`.

What COLD reproduces:

- **Exact** (`tests/test_cold_exact.py`, `tests/test_pfull_graph.py`, `tests/test_raw_graph.py`): the P_full
  graph rebuilt from the released inputs equals its recorded signature (`tests/pfull_rebuild_signature.json`:
  532,929 edges with the shipped broad graph, 528,008 with `DISATLAS_RAWGRAPH=1`; edge-set, row-order and `w_arm`
  SHA-256); the shipped broad graph (369,405 edges) and the signature of the released model's P_full (532,958
  edges, 13,042 nodes) are checked unchanged.
- **Parity** (`tests/test_cold_tolerance.py`): the retrained ensemble is compared with the released model by the
  pre-registered rule (5-seed mean endpoint cosine on the 1,094 IxIDN-Matched test pairs; 95 % disease-cluster
  bootstrap CI of the AUROC difference, B = 2,000; PARITY if the lower bound is > -0.02 and coverage >= 95 %).
  The authors' 5-seed retrain on P_full gives test AUROC 0.8589 against 0.8605 for the released model,
  ΔAUROC -0.002 [-0.008, +0.005], coverage 100 %: PARITY. `python scripts/pfull_parity.py` replays this from
  the shipped cosines (`tests/test_pfull_parity.py`). Walk sampling and floating-point sums depend on the
  thread count and BLAS build, so a retrain is judged by this rule, not by identity.
- The warm test suite then runs as in a warm run, so the manuscript numbers stay exact.

Measured on fresh clones, one cluster node each (peak memory is the Slurm MaxRSS of the job):

| Mode | Wall time | CPUs (memory allocated) | Peak memory |
|---|---|---|---|
| WARM | WARM_TIME_TBD | 16 (64 GB) | WARM_MEM_TBD |

## Data and licensing

Of the 918 curated disease–gene pairs used for the graph of record, the 604 supported by OMIM alone are not redistributed for licence reasons, and 6 more with no public source are not released. The 308 released pairs (`data/raw/public_disease_gene/`) are all in the public Open Targets 26.09 direct evidence, below the thresholds used to build the shipped Open Targets files. The remaining released BIKG-derived data are shared with permission from the AstraZeneca data office. All reported results reproduce exactly from the shipped graph files and model outputs. Rebuilding the gene-anchored graph from the public inputs here changes 36,868 edges (6.9%) relative to the training graph of the released model; a 5-seed retrain on a rebuild that differs from this one by 314 edges (0.06%) is within retraining noise (IxIDN-Matched test ΔAUROC −0.0016, 95% CI −0.0074 to +0.0041). The note "rebuilt from shipped inputs" in the shipped `model_config.json` predates the exclusion of the OMIM-only pairs; the final graph reproduces from the shipped graph files.

- `data/raw` holds minimal tables: only the columns the code reads, and only the rows that pass the published
  thresholds (no pre-filter data).
- Content from licence-restricted sources is blinded, not dropped: restricted targets and modules are replaced by
  stable salted tokens (`T_<12 hex>`, `M_<12 hex>`), so every overlap and graph computation is preserved exactly but
  the token cannot be mapped back. `disatlas.blinding.is_blinded(x)` tells you whether an identifier is a token.
- Results that depend on licensed trial and approval records are shipped as precomputed tables in `data/derived`.
- LLM-derived inputs (IxIDN-Matched curation verdicts, comparator system scores) are the authors' stored outputs,
  shipped in `data/raw/benchmarks` and used as-is; the pipeline calls no language model. The prompts are given in the
  manuscript and Supplementary Information, and you can run them with your own LLM on the blinded, shuffled benchmark
  files listed in [data/README.md, 'LLM-derived inputs'](data/README.md#llm-derived-inputs) and compare the results
  with the stored outputs.
- Column dictionaries, blinding details and the Zenodo file list: [data/README.md](data/README.md).
  Data licences: [data/README.md, 'Licence and blinding'](data/README.md#licence-and-blinding).

## Tests and parity

```bash
.venv/bin/python -m pytest -q tests           # or ./run_pipeline.sh --only tests (figures must be built first)
```

| Test | Checks |
|---|---|
| `tests/test_data_model.py` | blinding tokens, benchmark counts, broad-graph build statistics, ensemble reproduces the shipped IxIDN-Matched scores |
| `tests/test_pfull_graph.py` | P_full training graph rebuilt from the repository equals the recorded graph (edge set, row order, `w_arm`) |
| `tests/test_raw_graph.py` | broad graph rebuilt from `data/raw` equals the shipped one (`DISATLAS_RAWGRAPH=1`, about 1 h) |
| `tests/test_pfull_parity.py` | the recorded P_full retrain parity result replays from the shipped cosines |
| `tests/test_main_parity.py` | every number of the main text, from `tests/main_text_numbers.tsv` and the figure/table value files |
| `tests/test_sibench_parity.py` | numbers quoted in the SI benchmark sections |
| `tests/test_sibio_parity.py` | numbers quoted in the SI biology sections (specs in `tests/parity/`) |
| `tests/test_readme_links.py` | every relative link in this README (including the figure and table links) resolves to a tracked file |
| `tests/test_cold_exact.py` | COLD: the P_full graph in `build/graph` equals the recorded one |
| `tests/test_cold_tolerance.py` | COLD: the retrained ensemble passes the parity test |

A WARM run ends with `WARM_TESTS_TBD` and `DisAtlas pipeline finished (warm)`.

## Troubleshooting

- `data/zenodo missing: download the Zenodo bundle`: unpack the Zenodo record into `data/zenodo` (see Quick start).
- Checksum failure in `verify`: re-download the listed file; do not edit `data/ZENODO_SHA256SUMS.txt`.
- Memory: the all-pairs score matrix (13,021 x 13,021) is expanded in memory by several scripts; WARM peaks at about
  13 GB (measured 12.3-12.6 GiB in fresh-clone runs). On smaller machines run stages one at a time with `--only`.
- Threads: the retrain stage sets `OMP_NUM_THREADS` to the number of cores divided by five; lower it on shared machines.
- A figure script fails with a display error: set `MPLBACKEND=Agg` (the pipeline does this).

## How to cite

```bibtex
@article{disatlas,
  title   = {A Mechanistic Disease Representation Framework for Therapeutic Knowledge Transfer Across Rare and Common Conditions},
  author  = {AUTHORS_TBD},
  journal = {JOURNAL_TBD},
  year    = {YEAR_TBD},
  doi     = {DOI_TBD}
}
```

Data and model: Zenodo, DOI [10.5281/zenodo.24263240](https://doi.org/10.5281/zenodo.24263240).

## Licence

Code: MIT, see [LICENSE](LICENSE). Data: see [data/README.md, 'Licence and blinding'](data/README.md#licence-and-blinding).
