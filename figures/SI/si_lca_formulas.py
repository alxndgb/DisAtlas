"""Supplementary Figure SI_MechDist: ontology distance formulas (LCA distance across HPO and GO BP/MF/CC).

Appendix section: Threshold and Distance-Metric Sensitivity (#sec:SI_threshold_sensitivity), Figure #fig:SI_MechDist.
Formula sheet, no data inputs. d_max: HPO 25, GO BP 23, MF 18, CC 18; GO namespace weights 1.0 / 0.25 / 0.10;
GO hard hub-term mask at hub fraction 0.02 (pipeline as used). Optional IC term/gene weights and the soft hub penalty
are drawn for completeness and are not used in the pipeline (as stated in the caption).
Output: figures/SI/figS4_lca_similarity_diameter_revised.{png,pdf}
Run:    python figures/SI/si_lca_formulas.py
"""
import os, sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _style import save  # noqa: E402

RED = '#b03a48'     # annotation colour
BLUE = '#1f5fa8'    # optional weights and their effect
F = 17              # formula font size
N = 11.5            # note font size

# (y, kind, text) rows; kinds: 'f' formula, 'n' red note, 'new' blue note, 'sec' section label, 'rule'
rows = []
def f(t): rows.append(('f', t))
def n(t): rows.append(('n', t))
def new(t): rows.append(('new', t))
def sec(t): rows.append(('sec', t))
def rule(): rows.append(('rule', None))
def gap(h=0.5): rows.append(('gap', h))

sec('Ontology distance')
f(r'$Dist^{Onto}_{O}(t_1,t_2) = dist(t_1, lca_{t_1,t_2}) + dist(t_2, lca_{t_1,t_2})$')
n(r'$O \in \{HPO, BP, MF, CC\}$; $dist$ = number of is_a edges to the lowest common ancestor; precomputed LCA matrices, not changed by weighting')
rule()

sec('Optional term weights (shown for completeness; not used in the pipeline)')
f(r'$IC(t) = -\log\,\frac{n_t}{N}$,     $\widehat{IC}(t) = \max\left(\varepsilon,\ \frac{IC(t)}{\max_{t^\prime} IC(t^\prime)}\right)$')
n(r'$n_t$ = DGS diseases annotated with $t$ or any descendant of $t$ (true-path rule); $N$ = DGS diseases with $\geq 1$ annotation on that axis; floor $\varepsilon$ keeps all-hub diseases defined')
f(r'$w_{d,p} = \widehat{IC}(p)^{\gamma}\cdot e_{d,p}$          $w_{d,g} = e_{d,g}\ \ [\cdot\,\widehat{IC}(g)^{\gamma}]$,     $\gamma \in \{1, 2\}$')
n(r'$\gamma$ sharpens hub suppression: $\gamma=1$ gives developmental delay 0.15 vs cherry red spot 0.66 (4x); $\gamma=2$ gives 0.02 vs 0.44 (20x). Fever: 0.32 / 0.10')
f(r'$e_{d,p} = 0.5^{\,hops}\ \ [\cdot\,f_{d,p}]$          $e_{d,g} = s(c_{d,g})\ \ [\cdot\,q_{d,g}]$')
n(r'Symptom: $hops = 0$ for a direct HPOA annotation, $\geq 1$ when inherited from a child disease; $f_{d,p}$ = HPOA frequency')
n(r'Gene: $c_{d,g}$ = strongest source class of the edge (germline curated > OT aggregate > literature > GWAS solo, somatic solo); $s(c)$ = class prior, e.g. 1.0 / 1.0 / 1.0 / 0.5 / 0.5')
n(r'$q_{d,g}$ = score percentile within its own source (scores are not comparable across sources); $[\ ]$ = optional factor')
new(r'All $w = 1$  $\Rightarrow$  $\sum_i w_i = m$  $\Rightarrow$  the unweighted formulas exactly')
rule()

sec('HPO distance')
f(r'$Dist_p(p_1,p_2) = \frac{Dist^{Onto}_{HPO}(p_1,p_2)}{d^{HPO}_{max}}$')
n(r'$d^{HPO}_{max} = 25$, diameter of HPO')
gap(0.3)
f(r'$Dist^{HPO}(d_1,d_2) = \frac{1}{2}\left(\frac{\sum_{i=1}^{m} w_{d_1,p_i}\ \min_{1\leq j\leq n} Dist_p(p_i,p_j)}{\sum_{i=1}^{m} w_{d_1,p_i}} + \frac{\sum_{j=1}^{n} w_{d_2,p_j}\ \min_{1\leq i\leq m} Dist_p(p_i,p_j)}{\sum_{j=1}^{n} w_{d_2,p_j}}\right)$')
n(r'$m$, $n$ = number of $d_1$, $d_2$ symptoms; the weight belongs to the symptom being averaged, the best-match $\min$ is unweighted')
new(r'With weights, $\frac{1}{m}\sum_i$ and $\frac{1}{n}\sum_j$ become weighted means; a shared generic symptom (e.g. developmental delay) contributes a 0 with small weight')
rule()

sec('GO distance -- one per namespace (BP, MF, CC)')
for X, dm, lab in (('BP', 23, 'biological process'), ('MF', 18, 'molecular function'), ('CC', 18, 'cellular component')):
    f(r'$Dist^{%s}_g(g_1,g_2) = \frac{\min_{1\leq i\leq n^{%s}_{g_1},\ 1\leq j\leq m^{%s}_{g_2}}\ Dist^{Onto}_{%s}(g_1^i, g_2^j)}{d^{%s}_{max}}$' % (X, X, X, X, X))
    n(r'$n^{%s}_{g_1}$, $m^{%s}_{g_2}$ = number of GO %s terms of $g_1$, $g_2$;  $d^{%s}_{max} = %d$, diameter of GO %s' % (X, X, X, X, dm, lab))
n(r"Used in the pipeline: $n^{X}_{g_1}, m^{X}_{g_2}$ exclude hub-masked terms -- any GO term annotated to $>$2% of its namespace's DGS genes")
n(r'(e.g. "regulation of transcription by RNA pol II", on 900-1,100 DGS genes) is dropped from the gene-gene min before the formula runs; a masked-out gene keeps all its terms')
n(r'This hard mask (hub fraction 0.02) is applied to BP, MF and CC alike')
n(r'and is what produced the MechDist network used throughout')
f(r'$Dist^{X}_{g,\kappa}(g_1,g_2) = \min_{i,j}\ \left[(1-\kappa)\,\frac{Dist^{Onto}_{X}(g_1^i, g_2^j)}{d^{X}_{max}} + \kappa\left(1 - \widehat{IC}_X(lca_{ij})^{\gamma}\right)\right]$,     $\kappa = 0$ gives the hard-masked formula above')
n(r'Alternative shown for completeness (not used): a soft, graded penalty in place of the hard mask: a gene pair whose best match runs through a generic-but-unmasked GO term')
n(r'is downweighted rather than scored as identical. $\widehat{IC}_X$ is computed over DGS genes (true-path); $lca_{ij}$ is read from the precomputed LCA matrix; for an exactly shared term $lca = $ the term itself')
new(r'Hub GENES (e.g. TP53, on 152 diseases) would be handled at disease level by $\widehat{IC}(g)^{\gamma}$ in $w_{d,g}$ above (not used)')
new(r'Hub GO TERMS (this box) are removed by the hard mask in the pipeline; $\kappa$ is a finer alternative, not used')
rule()
gap(0.3)
f(r'$Dist^{X}(d_1,d_2) = \frac{1}{2}\left(\frac{\sum_{i=1}^{m} w_{d_1,g_i}\ \min_{1\leq j\leq n} Dist^{X}_g(g_i,g_j)}{\sum_{i=1}^{m} w_{d_1,g_i}} + \frac{\sum_{j=1}^{n} w_{d_2,g_j}\ \min_{1\leq i\leq m} Dist^{X}_g(g_i,g_j)}{\sum_{j=1}^{n} w_{d_2,g_j}}\right)$,   $X \in \{BP, MF, CC\}$')
n(r'$m$, $n$ = number of $d_1$, $d_2$ genes annotated in namespace $X$; weights do not touch gene-gene $Dist^X_g$, so gene-gene distances are unchanged')
f(r'$Dist^{GO}(d_1,d_2) = \sum_{X\,\in\,present} \lambda_X\, Dist^{X}(d_1,d_2)$,     $\lambda_{BP}=1.0,\ \lambda_{MF}=0.25,\ \lambda_{CC}=0.1$')
n(r'$present$ = namespaces in which both diseases have annotated genes (a weighted sum, not renormalised); NaN if none')
new(r'Optional (not used): gene weights $w_{d,g}$ in the average')
n(r'Weights are normalised per disease, so $s(c)$ matters only for diseases that mix source classes: 359 of 9,609 DGS diseases with genes (curated + GWAS/somatic)')

H = {'f': 0.95, 'n': 0.42, 'new': 0.42, 'sec': 0.55, 'rule': 0.35}
height = sum(H.get(k, v if k == 'gap' else 0) for k, v in rows) + 0.6
fig = plt.figure(figsize=(17, height))
y = height - 0.35
for k, t in rows:
    if k == 'gap':
        y -= t
        continue
    y -= H[k] / 2
    yy = y / height
    if k == 'f':
        fig.text(0.05, yy, t, fontsize=F, va='center')
    elif k == 'n':
        fig.text(0.07, yy, t, fontsize=N, color=RED, va='center')
    elif k == 'new':
        fig.text(0.07, yy, t, fontsize=N, color=BLUE, va='center', fontweight='bold')
    elif k == 'sec':
        fig.text(0.02, yy, t, fontsize=13.5, va='center', fontweight='bold')
    elif k == 'rule':
        fig.add_artist(plt.Line2D([0.02, 0.98], [yy, yy], color='k', lw=2.5, ls=(0, (6, 4))))
    y -= H[k] / 2
fig.text(0.98, 0.006, 'Black: formula; red: definitions and the pipeline as used; blue: optional weights and their effect. '
         'The GO hard hub-term mask is part of the pipeline; the term and gene weights and the soft penalty are shown for completeness and are not used.',
         ha='right', va='bottom', fontsize=9, color='#555')
print(save(fig, 'figS4_lca_similarity_diameter_revised', dpi=150))
