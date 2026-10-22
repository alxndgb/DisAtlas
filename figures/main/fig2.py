"""Figure 2 (Methods overview): DisAtlas pipeline schematic. Counts and coverage percentages in the boxes are
the values reported in the manuscript (Results and Methods, channel coverage).

    python figures/main/fig2.py      -> figures/main/fig2_pipeline.{png,pdf}
"""
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
mpl.rcParams.update({"font.family":"DejaVu Sans","font.size":7,"svg.fonttype":"none","pdf.fonttype":42})
GREY=("#ececec","#6b6b6b"); YEL=("#fdf2c8","#c9930e"); BLUE=("#e1ebf7","#3a64a0"); WHITE=("#ffffff","#7a2b2b"); PINK=("#f7e9e9","#7a2b2b")
T, B = 7.6, 6.8
fig = plt.figure(figsize=(7.2, 6.3), dpi=300)
ax = fig.add_axes([0,0,1,1]); ax.set_xlim(0,100); ax.set_ylim(0,87.5); ax.axis("off")

def box(x,y,w,h,title,body="",col=GREY,tfs=T):
    fc,ec=col
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0,rounding_size=0.8",fc=fc,ec=ec,lw=1.0))
    if body:
        ax.text(x+w/2,y+h-1.3,title,ha="center",va="top",fontsize=tfs,fontweight="bold")
        ax.text(x+w/2,y+h-4.0,body,ha="center",va="top",fontsize=B,linespacing=1.3)
    else:
        ax.text(x+w/2,y+h/2,title,ha="center",va="center",fontsize=tfs,fontweight="bold",linespacing=1.3)
def arr(p,q,**k):
    ax.add_patch(FancyArrowPatch(p,q,arrowstyle="-|>",mutation_scale=8,lw=1.0,color="#333333",shrinkA=0,shrinkB=0,**k))
def line(xs,ys): ax.plot(xs,ys,color="#333333",lw=1.0,solid_capstyle="butt")

# House panel letters (figures/main/disatlas_style.py): bold lower-case, sized so that
# at the same print width they match Figure 3 (15 pt in a 15.28 in wide figure -> 7.2 in here).
FIG3_WIDTH_IN = 4583 / 300          # rendered width of fig3_disease_map.png
PANEL_LABEL_PT = 15 * 7.2 / FIG3_WIDTH_IN   # = 7.07 pt
def plabel(x, y, letter):
    ax.text(x, y, letter, fontsize=PANEL_LABEL_PT, fontweight="bold", va="baseline", ha="left")

# ---------- Panel A ----------
plabel(0.8,84.55,"a")
ax.text(3.0,84.55,"Four evidence channels are fused into one disease graph and embedded",fontsize=8.5,va="baseline")
box(0.8,66.5,15,8,"Disease genes")
box(0.8,50.5,15,8,"Disease\nsymptoms")
ch=[("GO gene function","shared biological processes\n82% of pairs",GREY,74.0),
    ("HPA cell-type expression","shared cell-type profiles\n91% of pairs",YEL,65.0),
    ("Pfam protein domains","shared protein domains\n94% of pairs",YEL,56.0),
    ("HPO phenotype","shared symptoms\n35% of pairs",GREY,47.0)]
cx,cw,chh=21,24,8
for t,b,c,y in ch: box(cx,y,cw,chh,t,b,c)
for y in (74.0,65.0,56.0): arr((15.8,70.5),(cx,y+chh/2))
arr((15.8,54.5),(cx,47.0+chh/2))
fx,fy,fw,fh=50,57.5,19,15
box(fx,fy,fw,fh,"Rank-mean fusion","each pair is ranked\nin every channel;\nranks are averaged\n\nmissing channel\nranks last",YEL)
for t,b,c,y in ch: arr((cx+cw,y+chh/2),(fx,fy+fh/2))
rx,rw,rh=74,25,8.6
box(rx,72.0,rw,rh,"MechDist disease graph","13,021 diseases\nbenchmark pairs removed",BLUE)
box(rx,60.5,rw,rh,"node2vec embedding","32 dimensions\n= DisAtlas embeddings",BLUE)
box(rx,49.0,rw,rh,"Pair decoder","MLP, 5-seed ensemble\nscore for any disease pair",BLUE)
arr((fx+fw,fy+fh/2+3.0),(rx,76.3))
arr((rx+rw/2,72.0),(rx+rw/2,60.5+rh)); arr((rx+rw/2,60.5),(rx+rw/2,49.0+rh))
ax.text(cx+cw/2,45.6,"Missing annotations are not imputed",ha="center",va="top",fontsize=B,style="italic",color="#444444")

# ---------- Panel B ----------
plabel(0.8,39.95,"b")
ax.text(3.0,39.95,"An independent benchmark built from clinical trials",fontsize=8.5,va="baseline")
ax.add_patch(FancyBboxPatch((0.8,1.0),58.5,34.0,boxstyle="round,pad=0,rounding_size=1.2",fc="#f4f4f4",ec="#9a9a9a",lw=0.8))
ax.text(30,33.6,"Benchmark pairs from trial co-investigation",ha="center",va="top",fontsize=T,fontweight="bold",color="#333333")
box(2.5,11.5,13.5,12,"Drug–disease\ntrials","",GREY)
box(19.5,11.5,13.5,12,"Trial graph","diseases linked\nby shared drugs",YEL)
box(36.5,19.5,21.5,10.0,"Positive pairs (IxIDN⁺)","trialled with\nthe same drug",WHITE)
box(36.5,4.0,21.5,12.0,"Negative pairs (IxIDN⁻)","≥ 3 links apart, no\nshared-mechanism drug;\nmatched to positives",PINK)
arr((16,17.5),(19.5,17.5))
line([33,34.8],[17.5,17.5]); line([34.8,34.8],[10.0,24.5])
arr((34.8,24.5),(36.5,24.5)); arr((34.8,10.0),(36.5,10.0))
# benchmark + tasks
bx,bw=63,36
box(bx,5.5,bw,11.5,"IxIDN-Matched benchmark","disease-disjoint tune / test split\nmetric: AUROC",WHITE)
line([58,60.5],[24.5,24.5]); line([58,60.5],[10.0,10.0]); line([60.5,60.5],[10.0,24.5])
arr((60.5,11.25),(bx,11.25))
tw,th,ty=17,9.5,23.0
box(bx,ty,tw,th,"Task 1\nDisease clustering","",GREY,tfs=B+0.4)
box(bx+bw-tw,ty,tw,th,"Task 2\nDisease–disease\nassociation (DDAL)","",GREY,tfs=B+0.4)
arr((bx+tw/2,17.0),(bx+tw/2,ty)); arr((bx+bw-tw/2,17.0),(bx+bw-tw/2,ty))
ax.text(bx+bw/2,20.0,"evaluates",ha="center",va="center",fontsize=B,style="italic",color="#444444")
# embeddings feed tasks
ex=rx+rw/2
line([ex,ex],[49.0,36.5]); line([bx+tw/2,bx+bw-tw/2],[36.5,36.5])
arr((bx+tw/2,36.5),(bx+tw/2,ty+th)); arr((bx+bw-tw/2,36.5),(bx+bw-tw/2,ty+th))
ax.text(ex+0.8,42.5,"DisAtlas\nembeddings",ha="left",va="center",fontsize=B,style="italic",color="#444444")


if __name__ == "__main__":
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from disatlas_style import save
    save(fig, "fig2_pipeline", facecolor="white")
