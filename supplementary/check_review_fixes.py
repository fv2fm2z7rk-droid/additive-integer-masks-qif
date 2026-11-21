"""Exact-arithmetic checks of Lemma 15 (CSF version) and of the old greedy peeling.
(a) Lemma 15(1): V_p(T_{L,k}) = p + (1-p) g_k a^{k-1} = p + (2p-1) u_k, excess only at y = L+k.
(b) Lemma 15(2): Delta(-m) = g_k/(1-a) (2a^m - 1 - a^k); L=-m_k minimises the cost; C_k closed form;
    tie <=> L=-m_k-1 is also a minimiser.
(c) Historical: greedy peeling with alpha_i = q(L_i)(1-a^{k_i})/(1-a): reconstruction, sum alpha = 1, additivity of V_p and cost.
    The greedy peeling is NOT additive for V_p in general (see check_peeling_gap.py); the paper
    therefore uses the level-set decomposition of Lemma A, checked in levelset_check.py.
The definitions T, Vp, cost, mk, A, Ck are reused by verify_numbers.py.
"""
from fractions import Fraction as F
import random

def T(L,k,a):
    g=(1-a)/(1-a**k)
    return {L+j: g*a**j for j in range(k)}

def Vp(q,p):
    ys=set(q)|{y+1 for y in q}
    return p+sum(max(0,(1-p)*q.get(y-1,0)-p*q.get(y,0)) for y in ys)

def cost(q): return sum(abs(n)*w for n,w in q.items())

def mk(a,k): return sum(1 for j in range(1,k) if 2*a**j>1+a**k)
def A(a,m): return sum(a**i for i in range(1,m+1))
def Ck(a,k):
    m=mk(a,k); u=a**k/(1-a**k)
    return a/(1-a)+(m-2*A(a,m))+u*(2*(m-A(a,m))-k)

def peel(q,a):
    q=dict(q); comps=[]
    while any(v>0 for v in q.values()):
        L=min(n for n,v in q.items() if v>0); qL=q[L]; k=1
        while q.get(L+k,0)>=a**k*qL and q.get(L+k,0)>0: k+=1
        al=qL*(1-a**k)/(1-a)
        t=T(L,k,a)
        for n,w in t.items(): q[n]=q.get(n,0)-al*w
        assert all(v>=0 for v in q.values())
        comps.append((al,L,k))
    return comps

def main():
    from fractions import Fraction as F
    import random
    ok=True
    for p in [F(3,5),F(2,3),F(7,10),F(3,4),F(4,5),F(9,10),F(51,100),F(55,100)]:
        a=(1-p)/p
        for k in range(1,15):
            g=(1-a)/(1-a**k); u=a**k/(1-a**k)
            for L in range(-6,4):
                t=T(L,k,a)
                ok&= Vp(t,p)==p+(1-p)*g*a**(k-1)==p+(2*p-1)*u
            # (b)
            costs={L:cost(T(L,k,a)) for L in range(-k-3,4)}
            for m in range(0,k+1):
                ok&= costs[-m+1]-costs[-m]==g/(1-a)*(2*a**m-1-a**k)
            m=mk(a,k); cmin=min(costs.values())
            ok&= costs[-m]==cmin==Ck(a,k)
            tie=(2*a**(m+1)==1+a**k)
            ok&= ((costs[-m-1]==cmin)==tie)
    print("Lemma 15 (1),(2) checks:",ok)

    # tie example: a^{m+1} relation; p with a=1/2 -> k=1? search ties
    ties=[]
    for p in [F(2,3),F(3,5),F(5,8)]:
        a=(1-p)/p
        for k in range(1,12):
            m=mk(a,k)
            if 2*a**(m+1)==1+a**k: ties.append((p,k,m))
    print("ties found:",ties)

    # (c) peeling on random finitely supported masks
    random.seed(1)
    recon=sumone=addV=addC=True; nadd=0; N=300
    for _ in range(N):
        p=random.choice([F(3,5),F(7,10),F(4,5)]); a=(1-p)/p
        sup=range(random.randint(-4,0),random.randint(1,5))
        w=[F(random.randint(0,9)) for _ in sup]
        if sum(w)==0: continue
        q={n:x/sum(w) for n,x in zip(sup,w)}
        comps=peel(q,a)
        s=sum(al for al,_,_ in comps); sumone&=(s==1)
        V=sum(al*Vp(T(L,k,a),p) for al,L,k in comps)
        C=sum(al*cost(T(L,k,a)) for al,L,k in comps)
        addC&=(C==cost(q))
        if V==Vp(q,p): nadd+=1
    print("sum alpha == 1 always:",sumone," cost additive:",addC," V_p additive in",nadd,"of",N,"random masks")

if __name__ == "__main__":
    main()
