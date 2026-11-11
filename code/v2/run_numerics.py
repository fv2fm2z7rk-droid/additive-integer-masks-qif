"""Step 1 numerics: modular LP lower bounds (Wrapping lemma), rigorous integer
dual lower bounds, rigorous integer primal upper bounds, and the closed-form
envelope H_p of Theorem 2, on a grid of (p,c) with c in [c_lo(p), c_crit(p)).
c_lo = c_2(p) for p <= phi and 1-p for p > phi."""
import json, math, time, sys
import numpy as np
from fractions import Fraction as F
from core import int_primal, int_dual, rig_upper, rig_lower, mod_lp, c1c2, PHI
from frontier_exact import envelope_points, H, ccrit

PS = [F(52,100), F(53,100), F(55,100), F(57,100), F(3,5), F(618034,1000000), F(13,20),
      F(2,3), F(7,10), F(3,4), F(4,5), F(17,20), F(9,10), F(19,20)]
NS = [8, 16, 32, 64, 128, 256, 512]
NC = int(sys.argv[1]) if len(sys.argv) > 1 else 24
rows = []
t0 = time.time()
for pF in PS:
    p = float(pF); a = (1 - p) / p
    cc = ccrit(pF); ccf = float(cc)
    lo = c1c2(p)[1] if p <= PHI else 1 - p
    hull = envelope_points(pF, 400 if p < 0.56 else 200)
    lstar = math.ceil(math.log(2) / math.log(1 / a)) - 1
    A = 60
    B = int(min(900, max(200, lstar + 60 / (-math.log(a)))))
    for i in range(NC):
        c = lo + (ccf - lo) * (i + 0.5) / NC
        cF = F(c).limit_denominator(10**6); c = float(cF)
        vP, q, ys = int_primal(p, c, A, B)
        vD, al, lam, mu = int_dual(p, c, A, B)
        U, qx = rig_upper(pF, cF, q, ys)
        L, *_ = rig_lower(pF, cF, al, lam, A, B)
        h = H(pF, cF, hull=hull)
        mods = {n: mod_lp(p, c, n) for n in NS}
        nb = max(mods, key=mods.get)
        supp = [int(y) for y, v in zip(ys, q) if v > 1e-9]
        rows.append(dict(p=p, p_exact=str(pF), c=c, c_exact=str(cF), c_crit=ccf, A=A, B=B,
                         int_upper_rig=float(U), int_lower_rig=float(L), gap_rig=float(U - L),
                         H=float(h), U_minus_H=float(U - h), H_minus_L=float(h - L),
                         modular={str(n): v for n, v in mods.items()}, best_n=nb,
                         modular_best=mods[nb], U_minus_modular=float(U) - mods[nb],
                         support=[min(supp), max(supp)]))
    print(f"p={p:.4f} done ({time.time()-t0:.0f}s)", flush=True)
json.dump(dict(description=__doc__, rows=rows), open("data/numerics_grid.json", "w"), indent=1)
g = max(r['gap_rig'] for r in rows); gh = max(abs(r['U_minus_H']) for r in rows); gl = max(abs(r['H_minus_L']) for r in rows)
print("max rigorous gap U-L", g, " max |U-H|", gh, " max |H-L|", gl)
for pF in PS:
    rr = [r for r in rows if r['p_exact'] == str(pF)]
    print(f"p={float(pF):.4f}: max gap {max(r['gap_rig'] for r in rr):.2e}; max U - modular_best {max(r['U_minus_modular'] for r in rr):.2e}; worst c for modular {max(rr,key=lambda r:r['U_minus_modular'])['c']:.4f}")
