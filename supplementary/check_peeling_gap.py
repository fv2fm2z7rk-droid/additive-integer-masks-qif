"""Diagnose the peeling additivity claim: compare V_p(q) with sum_i alpha_i V_p(T_i),
and check that the frontier lower bound V_p(q) >= g(E|N|) still holds (Theorem 13 statement)."""
from fractions import Fraction as F
import random
from check_review_fixes import T,Vp,cost,peel,Ck
def g(c,p,K=60):
    a=(1-p)/p
    pts=[(Ck(a,k),p+(2*p-1)*a**k/(1-a**k)) for k in range(1,K)]
    for (c0,v0),(c1,v1) in zip(pts,pts[1:]):
        if c0<=c<=c1: return v0+(v1-v0)*(c-c0)/(c1-c0)
    return p
random.seed(1); bad=0; viol=0; ex=None
for _ in range(300):
    p=random.choice([F(3,5),F(7,10),F(4,5)]); a=(1-p)/p
    sup=range(random.randint(-4,0),random.randint(1,5))
    w=[F(random.randint(0,9)) for _ in sup]
    if sum(w)==0: continue
    q={n:x/sum(w) for n,x in zip(sup,w)}
    comps=peel(q,a)
    S=sum(al*Vp(T(L,k,a),p) for al,L,k in comps)
    if Vp(q,p)!=S:
        bad+=1
        if ex is None: ex=(p,q,comps,Vp(q,p),S)
    if Vp(q,p)<g(cost(q),p): viol+=1
print("non-additive:",bad,"  frontier violations:",viol)
p,q,comps,V,S=ex
print("example p=",p,"q=",{n:str(v) for n,v in q.items()})
print("components (alpha,L,k):",[(str(al),L,k) for al,L,k in comps])
print("V_p(q)=",V,float(V)," sum alpha V(T)=",S,float(S))
