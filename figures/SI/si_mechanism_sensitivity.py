
"""Supplementary Figure (fig:SI_mechanism_sensitivity): broad and specific evidence channels of the MechDist network.

  a  broad channel edges by evidence class    b  specific channel disease-module links by source
  c  combined graph edges by channel          d  specific edges by shared module source
  e  removing one channel or source and retraining (five seeds): test-split AUROC change and share of the
     50 nearest neighbours kept.
Bold panel letters aligned per column with headings beside them. Inputs: data/derived/mechanism_channels/
(channel_stats.json, broad_graph_size.json, broad_edges_by_class.tsv) and data/derived/ablation/.
Usage: python figures/SI/si_mechanism_sensitivity.py  (writes si_mechanism_sensitivity.png/.pdf next to this script)
"""
import sys as _s, pathlib as _p; _s.path.insert(0, str(_p.Path(__file__).resolve().parents[2] / "src")); import disatlas.deterministic  # noqa: E402,F401  reproducible PDF/SVG bytes
import os, pickle, json, collections, hashlib, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent; REPO = HERE.parents[1]
import sys as _sys; _sys.path.insert(0, str(REPO / "src"))
DV = REPO / "data" / "derived" / "mechanism_channels"; ABLD = REPO / "data" / "derived" / "ablation"
ABL = ABLD / "ablation_paired_delta_vs_released_B2000.tsv"; KNN = ABLD / "knn50_overlap_summary.tsv"
st = json.load(open(DV / "channel_stats.json"))
bg = json.load(open(DV / "broad_graph_size.json"))
# broad channel: edges by evidence class, counted on the built broad channel (before benchmark-pair stripping)
bc = pd.read_csv(DV / "broad_edges_by_class.tsv", sep="\t", index_col=0)["edges"].sort_values(ascending=False)
plt.rcParams.update({'font.size':7.5,'axes.spines.top':False,'axes.spines.right':False,'font.family':'DejaVu Sans',
                     'axes.labelsize':7.5,'xtick.labelsize':7,'ytick.labelsize':7,'axes.linewidth':0.7})
CB,CS,CR,CO='#4C72B0','#DD8452','#C44E52','#8172B3'
RL='Reactome pathway enrichment'
def nm(i): return RL if 'Reactome' in i else i.replace('Literature ','Lit. ')
fig=plt.figure(figsize=(7.2,8.4))
gs=fig.add_gridspec(3,2,height_ratios=[1,1,1.3],hspace=0.62,wspace=0.95,left=0.2,right=0.97,top=0.94,bottom=0.06)
P={}
# (a) broad
ax=fig.add_subplot(gs[0,0]); P['a']=(ax,'Broad channel: edges by evidence class',f"{bg['broad_graph_nodes']:,} diseases, {bg['broad_graph_edges']:,} edges")
LB={'hpo+px':'HPO + HPA/Pfam','go+px':'GO + HPA/Pfam','hpo_only':'HPO only','all':'GO + HPO + HPA/Pfam','go_only':'GO only','go+hpo':'GO + HPO'}
ax.barh(range(len(bc)),bc.values/1e3,color=CB,height=0.7); ax.set_yticks(range(len(bc))); ax.set_yticklabels([LB.get(i,str(i)) for i in bc.index]); ax.invert_yaxis()
for i,v in enumerate(bc.values): ax.text(v/1e3+3,i,f'{v:,}',va='center',fontsize=6.5)
ax.set_xlim(0,bc.max()/1e3*1.28); ax.set_xlabel('edges (thousands)')
# (b) specific links
t=pd.DataFrame(st['spec_by_type']).T.sort_values('links',ascending=False)
ax=fig.add_subplot(gs[0,1]); P['b']=(ax,'Specific channel: disease–module links by source',f"{st['spec_link_diseases']:,} diseases, {st['spec_modules']:,} modules")
ax.barh(range(len(t)),t.links/1e3,color=[CR if 'Reactome' in i else CS for i in t.index],height=0.7); ax.set_yticks(range(len(t))); ax.set_yticklabels([nm(i) for i in t.index]); ax.invert_yaxis()
for i,(l,m) in enumerate(zip(t.links,t.modules)): ax.text(l/1e3+3,i,f'{int(m):,} modules',va='center',fontsize=6.5)
ax.set_xlim(0,t.links.max()/1e3*1.45); ax.set_xlabel('disease–module links (thousands)')
# (c) overlap
ax=fig.add_subplot(gs[1,0]); P['c']=(ax,'Combined graph: edges by channel',f"{st['union_edges']:,} edges on {st['union_nodes']:,} diseases; {st['spec_nodes']:,} carry specific edges")
vals=[st['broad_only_edges'],st['overlap_edges'],st['spec_only_edges']]
ax.bar(['broad only','both','specific only'],np.array(vals)/1e3,color=[CB,CO,CS],width=0.6)
for i,v in enumerate(vals): ax.text(i,v/1e3+6,f'{v:,}',ha='center',fontsize=6.5)
ax.set_ylim(0,max(vals)/1e3*1.12); ax.set_ylabel('edges (thousands)')
# (d) shared source
s=pd.Series(st['spec_edges_by_shared_type']).sort_values(ascending=False)
ax=fig.add_subplot(gs[1,1]); P['d']=(ax,'Specific edges by shared module source','one edge can count for several sources')
ax.barh(range(len(s)),s.values/1e3,color=[CR if 'Reactome' in i else CS for i in s.index],height=0.7); ax.set_yticks(range(len(s))); ax.set_yticklabels([nm(i) for i in s.index]); ax.invert_yaxis()
for i,v in enumerate(s.values): ax.text(v/1e3+3,i,f'{int(v):,}',va='center',fontsize=6.5)
ax.set_xlim(0,s.max()/1e3*1.25); ax.set_xlabel('specific edges (thousands)')
# (e) removals
A=pd.read_csv(ABL,sep='\t'); A=A[A.case_set=='test'].set_index('arm'); K=pd.read_csv(KNN,sep='\t').set_index('arm')
arms=[('broad_off','Broad channel',CB,True),('broad_noGO','GO',CB,False),('broad_noHPO','HPO',CB,False),('broad_noHPA','HPA',CB,False),('broad_noPfam','Pfam',CB,False),
      ('spec_off','Specific channel',CS,True),('drop_REACTOME','Reactome modules',CR,False),('drop_LITGO','Lit. GO processes',CS,False),('drop_LITMECH','Lit. mechanism entities',CS,False),
      ('drop_LITPW','Lit. pathways',CS,False),('drop_NAMED','Curated gene sets',CS,False),('drop_lesion','Pathology text',CS,False)]
gs2=gs[2,:].subgridspec(1,2,wspace=0.08,width_ratios=[1.15,1])
ax1=fig.add_subplot(gs2[0]); ax2=fig.add_subplot(gs2[1],sharey=ax1)
P['e']=(ax1,'Removing one channel or source and retraining (five seeds)',None)
for i,(a,l,c,_) in enumerate(arms):
    r=A.loc[a]; ax1.errorbar(r.d_auroc,i,xerr=[[r.d_auroc-r.d_auroc_ci_lo],[r.d_auroc_ci_hi-r.d_auroc]],fmt='o',color=c,ms=3.8,capsize=2,lw=1)
    ax2.plot(K.loc[a,'mean'],i,'o',color=c,ms=3.8)
for ax_ in (ax1,ax2): ax_.axhline(4.5,color='0.85',lw=0.6)
ax1.axvline(0,color='grey',lw=0.6,ls='--'); ax1.set_yticks(range(len(arms))); ax1.set_yticklabels([l if g else '   '+l for _,l,_,g in arms]); ax1.invert_yaxis()
for tl,(_,_,_,g) in zip(ax1.get_yticklabels(),arms):
    if g: tl.set_fontweight('bold')
ax1.set_xlabel('change in test-split AUROC (95% CI)'); ax1.set_xticks([-0.1,-0.05,0])
sb=K.loc['A_seed_baseline','mean']; ax2.axvline(sb,color='k',lw=0.8); ax2.axvline(50/13020,color='grey',lw=0.6,ls=':')
ax2.text(sb+0.015,1.5,f'two seeds\nof released\nmodel ({sb:.2f})',fontsize=6.3,ha='left',va='center')
ax2.set_xlim(0,0.95); plt.setp(ax2.get_yticklabels(),visible=False); ax2.tick_params(axis='y',length=0)
ax2.set_xlabel('share of 50 nearest neighbours kept')
ax1.set_title('AUROC change',loc='center',fontsize=7,color='0.3'); ax2.set_title('Neighbourhood kept',loc='center',fontsize=7,color='0.3')
# panel letters: one x per column (left edge of the widest tick labels in that column), heading beside the letter
fig.canvas.draw(); REN=fig.canvas.get_renderer(); inv=fig.transFigure.inverted()
def tb(ax): return inv.transform(ax.get_tightbbox(REN))
colx={0:min(tb(P[k][0])[0,0] for k in 'ace'), 1:min(tb(P[k][0])[0,0] for k in 'bd')}
colx[0]=min(colx[0],0.012)
for k,(ax,h,sub) in P.items():
    x=colx[0] if k in 'ace' else colx[1]
    y=ax.get_position().y1+(0.046 if sub else 0.03)
    if k=='e': y=max(ax1.get_position().y1,ax2.get_position().y1)+0.032
    fig.text(x,y,k,fontsize=10,fontweight='bold',va='baseline',ha='left')
    fig.text(x+0.03,y,h,fontsize=8,fontweight='bold',va='baseline',ha='left')
    if sub: fig.text(x+0.03,y-0.018,sub,fontsize=7,color='0.3',va='baseline',ha='left')
fig.savefig(HERE / 'si_mechanism_sensitivity.png', dpi=300); fig.savefig(HERE / 'si_mechanism_sensitivity.pdf')
