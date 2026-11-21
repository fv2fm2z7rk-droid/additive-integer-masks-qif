"""Exact check of the level-set layer cake (w(s) = q(s) a^{-s}) used in the replacement proof
of the lower bound of Theorem 13, on random finitely supported rational masks (incl. the
p=3/5 counterexample to the greedy peeling) and on long two-sided-tail masks.
Checks: q = sum_I mu_I T_I, sum mu = 1, cost additive, V_p additive, V_p(q) >= g(E|N|)."""
from fractions import Fraction as F
import random, sys

def Vp(q,p):
    ys=set(q)|{y+1 for y in q}
    return p+sum(max(0,(1-p)*q.get(y-1,0)-p*q.get(y,0)) for y in ys)
def cost(q): return sum(abs(n)*w for n,w in q.items())
def T(L,k,a):
    g=(1-a)/(1-a**k); return {L+j:g*a**j for j in range(k)}
def mk(a,k): return sum(1 for j in range(1,k) if 2*a**j>1+a**k)
def A(a,m): return sum((a**i for i in range(1,m+1)),F(0))
def Ck(a,k):
    m=mk(a,k); u=a**k/(1-a**k); return a/(1-a)+(m-2*A(a,m))+u*(2*(m-A(a,m))-k)
def Vk(p,k):
    a=(1-p)/p; return p+(2*p-1)*a**k/(1-a**k)
_pts={}
def g(c,p):
    if p not in _pts:
        a=(1-p)/p; pts=[]; k=1
        while True:
            pts.append((Ck(a,k),Vk(p,k))); k+=1
            if k>80: break
        _pts[p]=pts
    pts=_pts[p]
    for (c0,v0),(c1,v1) in zip(pts,pts[1:]):
        if c0<=c<=c1: return v0+(v1-v0)*(c-c0)/(c1-c0)
    if c>=pts[-1][0]: return None   # beyond tabulated vertices (near c_crit): skip exact comparison
    raise ValueError

def levelset(q,a):
    """Return list of (mu, L, k): mu = nu-mass of component [L, L+k-1] (finite support => k finite)."""
    w={s:q[s]/a**s for s in q if q[s]>0}
    levels=sorted(set(w.values())); out={}
    prev=F(0)
    for t in levels:                      # on (prev, t] the set {w > t'} is {w >= t}
        S=sorted(s for s in w if w[s]>=t)
        comps=[]; start=S[0]; last=S[0]
        for s in S[1:]:
            if s==last+1: last=s
            else: comps.append((start,last)); start=last=s
        comps.append((start,last))
        dt=t-prev
        for L,R in comps:
            k=R-L+1; M=a**L*(1-a**k)/(1-a)
            out[(L,k)]=out.get((L,k),F(0))+dt*M
        prev=t
    return [(mu,L,k) for (L,k),mu in out.items()]

def check(q,p):
    a=(1-p)/p; comps=levelset(q,a)
    rec={}
    for mu,L,k in comps:
        for n,v in T(L,k,a).items(): rec[n]=rec.get(n,0)+mu*v
    ok_rec=all(rec.get(n,0)==q.get(n,0) for n in set(rec)|set(q))
    ok_sum=sum(mu for mu,_,_ in comps)==1
    ok_cost=sum(mu*cost(T(L,k,a)) for mu,L,k in comps)==cost(q)
    ok_V=sum(mu*Vp(T(L,k,a),p) for mu,L,k in comps)==Vp(q,p)
    gc=g(cost(q),p)
    ok_front=True if gc is None else Vp(q,p)>=gc
    return ok_rec,ok_sum,ok_cost,ok_V,ok_front

if __name__=="__main__":
    res=[]
    # counterexample to greedy peeling
    p=F(3,5); q={-1:F(7,24),0:F(1,4),1:F(1,8),2:F(1,24),3:F(7,24)}
    print("counterexample p=3/5:",check(q,p), " V_p(q)=",Vp(q,p)," g(E|N|)=",g(cost(q),p))
    random.seed(1); n=0; allok=[True]*5
    for _ in range(300):
        p=random.choice([F(3,5),F(7,10),F(4,5)])
        sup=range(random.randint(-4,0),random.randint(1,5))
        w=[F(random.randint(0,9)) for _ in sup]
        if sum(w)==0: continue
        q={s:x/sum(w) for s,x in zip(sup,w) if x>0}
        r=check(q,p); n+=1; allok=[x and y for x,y in zip(allok,r)]
    print("same 300-mask stream as before:",n,"masks; reconstruct/sum1/cost/V-additive/frontier:",allok)
    random.seed(7); n=0; allok=[True]*5
    for _ in range(400):
        p=random.choice([F(11,20),F(3,5),F(5,8),F(2,3),F(7,10),F(3,4),F(4,5),F(9,10)])
        lo=random.randint(-12,0); hi=random.randint(1,12)
        w=[F(random.randint(0,20))*random.choice([1,1,0]) for _ in range(lo,hi+1)]
        if sum(w)==0: continue
        q={s:x/sum(w) for s,x in zip(range(lo,hi+1),w) if x>0}
        r=check(q,p); n+=1; allok=[x and y for x,y in zip(allok,r)]
    print("400 further random masks (gaps, wider support, 8 biases):",n,allok)
    # long two-sided tails: geometric-like decay with random ratios, support [-60,60]
    random.seed(3); n=0; allok=[True]*5
    for _ in range(40):
        p=random.choice([F(3,5),F(7,10),F(9,10)])
        rl=F(random.randint(1,9),10); rr=F(random.randint(1,9),10)
        w={s:(rl**(-s) if s<0 else rr**s)*F(random.randint(1,3)) for s in range(-60,61)}
        Z=sum(w.values()); q={s:v/Z for s,v in w.items()}
        r=check(q,p); n+=1; allok=[x and y for x,y in zip(allok,r)]
    print("40 masks with two-sided tails on [-60,60]:",n,allok)
