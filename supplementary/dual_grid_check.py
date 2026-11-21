"""Rigorous LP checks of the Theorem 13 frontier on a grid of (p, c), reusing (read-only) the
certificate routines of core.py (in this folder):
 - lower: weak-duality certificate alpha:Z->[0,1] (alpha=1 left of window, 0 right), rationalised,
          mu recomputed exactly over all y in Z (tails included) => exact lower bound on V*_p(c)
          valid for masks of ANY support;
 - upper: window primal LP mask, rationalised => exact upper bound.
Compared with the closed form g(c) in exact arithmetic."""
from fractions import Fraction as F
import core
from levelset_check import g
ps=[F(11,20),F(3,5),F(618034,1000000),F(2,3),F(7,10),F(3,4),F(4,5),F(9,10)]
fr=[F(1,20),F(1,5),F(7,20),F(1,2),F(13,20),F(4,5),F(9,10),F(19,20)]
worst_lo=worst_up=F(0); viol=0; rows=[]
for p in ps:
    cc=F(core.c_crit_exact(p)).limit_denominator(10**15)
    W=40 if p>=F(3,5) else 70
    for f in fr:
        c=(f*cc).limit_denominator(10**6)
        gc=g(c,p)
        if gc is None: continue
        val,alpha,lam,mu=core.int_dual(float(p),float(c),W,W)
        lo,_,_,_=core.rig_lower(p,c,alpha,lam,W,W)
        up_f,qf,ys=core.int_primal(float(p),float(c),W,W)
        up,_=core.rig_upper(p,c,qf,ys)
        if lo>gc or up<gc: viol+=1
        worst_lo=max(worst_lo,gc-lo); worst_up=max(worst_up,up-gc)
        rows.append((float(p),float(c),float(lo),float(gc),float(up)))
for r in rows: print("p=%.6f c=%.6f  lower=%.12f  g=%.12f  upper=%.12f"%r)
print("points:",len(rows)," violations (lower>g or upper<g):",viol,
      " max g-lower=%.2e  max upper-g=%.2e"%(float(worst_lo),float(worst_up)))
