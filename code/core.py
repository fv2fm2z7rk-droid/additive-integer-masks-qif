"""Core LP / certificate routines for the biased-bit integer frontier V*_p(c).

V_p(q) = sum_s max{p q(s), (1-p) q(s-1)},  cost E|Y| <= c,  Y on Z.
Upper bounds: explicit masks (exactly evaluated in rationals).
Lower bounds: weak-duality certificates alpha: Z->[0,1] with alpha=1 left of the
window and 0 right of it (exactly checked in rationals, tails included), and
modular LP values on Z_n (Wrapping lemma), the latter in floating point.
"""
from fractions import Fraction as F
import math
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix, csr_matrix

PHI = (math.sqrt(5) - 1) / 2

def a_of(p): return (1 - p) / p

def c_crit(p):
    a = a_of(p)
    l = math.ceil(math.log(2) / math.log(1 / a)) - 1
    best = min(l0 + a * (2 * a ** l0 - 1) / (1 - a) for l0 in range(max(0, l - 2), l + 3))
    return best

def c_crit_exact(p):
    a = (1 - p) / p
    vals = []
    l = 0
    while True:
        vals.append(l + a * (2 * a ** l - 1) / (1 - a))
        if 2 * a ** (l + 1) < 1 and l > 0 and vals[-1] > vals[-2]:
            break
        l += 1
        if l > 10000: break
    return min(vals)

def c1c2(p):
    a = a_of(p)
    return (1 + a * a) / (1 + a + a * a), (1 + a * a + 2 * a ** 3) / (1 + a + a * a + a ** 3)

def known_V(p, c):
    """Integer paper I1-I3 closed forms where proven, else None."""
    if c <= 1 - p: return 1 - p * c
    if c >= c_crit(p): return p
    if p <= PHI:
        c1, c2 = c1c2(p)
        if c <= c1: return p + (1 - p) ** 2 / p * (1 - c)
        if c <= c2: return 1 - p + p ** 3 * (2 - c)
    return None

# ---------------------------------------------------------------- primal (upper bound)
def int_primal(p, c, A, B):
    """min V_p(q) over q on {-A..B}, E|Y|<=c.  Returns (fun, q, ys)."""
    ys = np.arange(-A, B + 1); N = len(ys)
    S = np.arange(-A, B + 2); M = len(S)
    nv = N + M
    rows = lil_matrix((2 * M + 1, nv)); b = np.zeros(2 * M + 1); r = 0
    for si, s in enumerate(S):
        # e_s >= p q(s) - and the sum of t_s = p + sum e_s with e_s >= (1-p)q(s-1)-p q(s)
        y = s - 1
        if -A <= y <= B:
            rows[r, y + A] = (1 - p)
            if s <= B: rows[r, s + A] = -p
            rows[r, N + si] = -1; r += 1
    rows[r, :N] = np.abs(ys); b[r] = c; r += 1
    Aub = csr_matrix(rows[:r]); b = b[:r]
    cobj = np.r_[np.zeros(N), np.ones(M)]
    Aeq = csr_matrix(np.r_[np.ones(N), np.zeros(M)].reshape(1, -1))
    res = linprog(cobj, A_ub=Aub, b_ub=b, A_eq=Aeq, b_eq=[1], bounds=[(0, None)] * nv,
                  method="highs", options=dict(primal_feasibility_tolerance=1e-10,
                                               dual_feasibility_tolerance=1e-10))
    assert res.status == 0, res.message
    return p + res.fun, res.x[:N], ys

def exact_V(pF, q):
    """q: dict y->Fraction (prob. law). exact V_p."""
    keys = set(q) | {y + 1 for y in q}
    tot = F(0)
    for s in keys:
        tot += max(pF * q.get(s, F(0)), (1 - pF) * q.get(s - 1, F(0)))
    return tot

def rig_upper(pF, cF, qfloat, ys, den=10**12):
    """Rationalise an LP mask, renormalise, and if the exact cost exceeds c
    mix with delta_0.  Returns an exact rational upper bound on V*_p(c)."""
    q = {}
    for y, v in zip(ys, qfloat):
        if v > 1e-15:
            q[int(y)] = F(int(round(v * den)), den)
    tot = sum(q.values()); q = {y: v / tot for y, v in q.items()}
    cost = sum(abs(y) * v for y, v in q.items())
    if cost > cF:
        th = 1 - cF / cost
        q = {y: (1 - th) * v for y, v in q.items()}
        q[0] = q.get(0, F(0)) + th
    assert sum(q.values()) == 1 and sum(abs(y) * v for y, v in q.items()) <= cF
    return exact_V(pF, q), q

# ---------------------------------------------------------------- dual (lower bound)
def int_dual(p, c, A, B):
    """max mu - lam c  s.t.  kappa(y) >= mu - lam|y| for all y in Z,
    alpha_y = 1 (y <= -A), alpha_y = 0 (y >= B+1), alpha free in [0,1] on [-A+1, B].
    Tails: y <= -A-1 -> kappa = p ; y >= B+1 -> kappa = 1-p.
    Variables: alpha_{-A+1..B} (N), lam, mu."""
    al_idx = {y: i for i, y in enumerate(range(-A + 1, B + 1))}
    N = len(al_idx); il, im = N, N + 1; nv = N + 2
    rows = lil_matrix((B + A + 3, nv)); b = []; r = 0
    def alpha_terms(y):  # returns (const, {idx: coeff}) for kappa(y)
        const = 0.0; co = {}
        # p alpha_y
        if y <= -A: const += p
        elif y <= B: co[al_idx[y]] = co.get(al_idx[y], 0) + p
        # (1-p)(1 - alpha_{y+1})
        z = y + 1
        if z <= -A: pass
        elif z <= B: const += (1 - p); co[al_idx[z]] = co.get(al_idx[z], 0) - (1 - p)
        else: const += (1 - p)
        return const, co
    for y in range(-A - 1, B + 2):
        const, co = alpha_terms(y)
        # mu - lam|y| - kappa <= 0  ->  mu - lam|y| - sum co*alpha <= const
        for i, v in co.items(): rows[r, i] = -v
        rows[r, il] = -abs(y); rows[r, im] = 1; b.append(const); r += 1
    Aub = csr_matrix(rows[:r])
    cobj = np.zeros(nv); cobj[im] = -1; cobj[il] = c
    bounds = [(0, 1)] * N + [(0, None), (None, None)]
    res = linprog(cobj, A_ub=Aub, b_ub=np.array(b), bounds=bounds, method="highs",
                  options=dict(primal_feasibility_tolerance=1e-10, dual_feasibility_tolerance=1e-10))
    assert res.status == 0, res.message
    alpha = {y: res.x[i] for y, i in al_idx.items()}
    return -res.fun, alpha, res.x[il], res.x[im]

def rig_lower(pF, cF, alpha, lam, A, B, den=10**12):
    """Exact lower bound from a rationalised certificate (all y in Z checked)."""
    al = {y: min(F(1), max(F(0), F(int(round(v * den)), den))) for y, v in alpha.items()}
    def A_(y):
        if y <= -A: return F(1)
        if y > B: return F(0)
        return al[y]
    lamF = max(F(0), F(int(round(lam * den)), den))
    # mu = inf_y kappa(y) + lam|y|; for y <= -A-1 it is p + lam|y| (min at y=-A-1),
    # for y >= B+1 it is 1-p+lam*y (min at y=B+1).
    vals = [pF * A_(y) + (1 - pF) * (1 - A_(y + 1)) + lamF * abs(y) for y in range(-A - 1, B + 2)]
    mu = min(vals)
    return mu - lamF * cF, al, lamF, mu

# ---------------------------------------------------------------- modular LP on Z_n
def mod_lp(p, c, n):
    """V*_{p,n}(c): min sum_s max{p q(s),(1-p) q(s-1)} on Z_n, E d(Y,0)<=c."""
    d = np.array([min(y, n - y) for y in range(n)], float)
    nv = 2 * n
    rows = lil_matrix((n + 1, nv)); b = np.zeros(n + 1)
    for s in range(n):
        rows[s, (s - 1) % n] += (1 - p); rows[s, s] += -p; rows[s, n + s] = -1
    rows[n, :n] = d; b[n] = c
    cobj = np.r_[np.zeros(n), np.ones(n)]
    Aeq = csr_matrix(np.r_[np.ones(n), np.zeros(n)].reshape(1, -1))
    res = linprog(cobj, A_ub=csr_matrix(rows), b_ub=b, A_eq=Aeq, b_eq=[1],
                  bounds=[(0, None)] * nv, method="highs",
                  options=dict(primal_feasibility_tolerance=1e-10, dual_feasibility_tolerance=1e-10))
    assert res.status == 0
    return p + res.fun
