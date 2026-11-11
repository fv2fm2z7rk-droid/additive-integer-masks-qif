"""Smallest modulus n for which the modular LP value V*_{p,n}(c) (a lower bound on
V*_p(c) by the Wrapping lemma) matches the closed-form integer frontier H_p(c)
to 1e-9.  Floating-point LP (HiGHS); numerical evidence only."""
import json
from fractions import Fraction as F
from core import mod_lp
from frontier_exact import envelope_points, H, Qk, ccrit, m_of
res = []
for pF in [F(11,20), F(3,5), F(618034,1000000), F(7,10), F(3,4), F(9,10)]:
    p = float(pF); hull = envelope_points(pF, 200); a = (1-pF)/pF
    for k in range(2, 13):
        Ck, Ck1 = Qk(pF, k)[0], Qk(pF, k+1)[0]
        c = (Ck + Ck1) / 2
        h = float(H(pF, c, hull=hull))
        nmin = None
        for n in range(2, 200):
            if h - p < 1e-6: break
            if mod_lp(p, float(c), n) >= h - 1e-9: nmin = n; break
        res.append(dict(p=p, k=k, c=float(c), H=h, n_min=nmin, m_k=m_of(a,k), m_k1=m_of(a,k+1)))
        print(f"p={p:.4f} cell k={k} c={float(c):.4f} n_min={nmin}  (m_k,m_k+1)=({m_of(a,k)},{m_of(a,k+1)})", flush=True)
json.dump(dict(description=__doc__, rows=res), open("data/modular_n_needed.json", "w"), indent=1)
