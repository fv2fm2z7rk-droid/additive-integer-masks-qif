"""Recompute every hand-typed number of the CSF paper from the definitions, in exact arithmetic.

Covers: Example 14 (p = 0.7), the appendix table (p = 0.6 and p = 0.8, k = 1..5), the masses in
the caption of Fig. 3 (p = phi, exact in Q(sqrt 5) via sympy), and the c_1, c_2 formulas and
explicit optimisers of Corollary 21, and the tail formulas of Corollary 20.  Every C_k and V_k is compared with a brute-force minimum of
the cost of T_{L,k} over all shifts L (and with V_p evaluated directly on the minimiser), not only
with the closed forms of Definition 12.  Also checks the increment formulas of Lemma 16, the
range of theta_k, the claim that k-1 and k are never both transitions, that every tie is a
transition, and the bound "at most l* ties" of Corollary 18.

Run from this folder:  python3 verify_numbers.py   (prints PASS/FAIL per item and overall)
"""
from fractions import Fraction as F
import math
import sympy as sp
from check_review_fixes import T, Vp, cost, mk, A, Ck

RESULTS = []
def check(name, cond):
    RESULTS.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL") + "  " + name)

def Vk_closed(p, k):
    a = (1 - p) / p
    return p + (2 * p - 1) * a**k / (1 - a**k)

def brute(p, k):
    """min over shifts L of the cost of T_{L,k}; returns (min cost, set of minimising L, V_p at minimisers)."""
    a = (1 - p) / p
    costs = {L: cost(T(L, k, a)) for L in range(-k - 3, 4)}
    cmin = min(costs.values())
    Ls = sorted(L for L, c in costs.items() if c == cmin)
    Vs = {Vp(T(L, k, a), p) for L in Ls}
    # minimiser must lie strictly inside the scanned window
    assert min(Ls) > -k - 3 and max(Ls) < 3
    return cmin, Ls, Vs

def lstar(a):
    j = 0
    while a**(j + 1) > F(1, 2): j += 1
    return j

def ccrit(p):
    a = (1 - p) / p
    return min(l + a * (2 * a**l - 1) / (1 - a) for l in range(0, lstar(a) + 3))

def theta(a, k):
    m, m1 = mk(a, k), mk(a, k + 1)
    if m1 == m: return 0 * a
    return (2 * a**(m + 1) - 1 - a**(k + 1)) / (a**k * (1 - a))

def Ek_formula(a, k):
    m = mk(a, k)
    return k - A(a, k) - 2 * (m - A(a, m)) - theta(a, k) * (1 - a**k)

def fmt4(x):
    return "%.4f" % float(x)

# ------------------------------------------------------------------ Example 14 (p = 0.7)
print("== Example 14, p = 7/10")
p = F(7, 10); a = (1 - p) / p
check("a = 3/7", a == F(3, 7))
check("l* = #{j>=1: a^j > 1/2} = 0", lstar(a) == 0)
check("c_crit(0.7) = a/(1-a) = 3/4", ccrit(p) == a / (1 - a) == F(3, 4))
expQ = {1: (F(0), F(1)), 2: (F(3, 10), F(79, 100)), 3: (F(39, 79), F(58, 79)), 4: (F(177, 290), F(4141, 5800))}
for k, (c, v) in expQ.items():
    cb, Ls, Vs = brute(p, k)
    check(f"Q_{k} = ({c}, {v}): closed form", (Ck(a, k), Vk_closed(p, k)) == (c, v))
    check(f"Q_{k}: brute force over shifts (minimising L = {Ls})", cb == c and Vs == {v})
check("Q_3 ~ (0.4937, 0.7342)", (fmt4(F(39, 79)), fmt4(F(58, 79))) == ("0.4937", "0.7342"))
check("Q_4 ~ (0.6103, 0.7140)", (fmt4(F(177, 290)), fmt4(F(4141, 5800))) == ("0.6103", "0.7140"))
check("m_k = 0 for k <= 40 at p = 0.7 (optimal shift L = 0)", all(mk(a, k) == 0 for k in range(1, 41)))

# ------------------------------------------------------------------ Appendix table
print("== Appendix table")
TABLE = {  # p: {k: (C_k exact, V_k exact, C_k shown, V_k shown)}
    F(3, 5): {1: (F(0), F(1), "0", "1"), 2: (F(2, 5), F(19, 25), "0.4000", "0.7600"),
              3: (F(13, 19), F(13, 19), "0.6842", "0.6842"), 4: (F(11, 13), F(211, 325), "0.8462", "0.6492"),
              5: (F(213, 211), F(133, 211), "1.0095", "0.6303")},
    F(4, 5): {1: (F(0), F(1), "0", "1"), 2: (F(1, 5), F(21, 25), "0.2000", "0.8400"),
              3: (F(2, 7), F(17, 21), "0.2857", "0.8095"), 4: (F(27, 85), F(341, 425), "0.3176", "0.8024"),
              5: (F(112, 341), F(273, 341), "0.3284", "0.8006")},
}
for p, rows in TABLE.items():
    a = (1 - p) / p
    for k, (c, v, cs, vs) in rows.items():
        cb, Ls, Vs = brute(p, k)
        ok = Ck(a, k) == c and Vk_closed(p, k) == v and cb == c and Vs == {v}
        ok &= (k == 1) or (fmt4(c) == cs and fmt4(v) == vs)
        check(f"p={p} k={k}: C={c} ({cs}), V={v} ({vs}) [closed form = brute force]", ok)
    check(f"p={p}: V_2 = 1-p+p^2 and C_2 = 1-p (Thm 10)", Vk_closed(p, 2) == 1 - p + p * p and Ck(a, 2) == 1 - p)
    Cs = [Ck(a, k) for k in range(1, 61)]; Vs_ = [Vk_closed(p, k) for k in range(1, 61)]
    cc = ccrit(p)
    check(f"p={p}: V_k strictly decreasing, > p, V_60 - p < 1e-6",
          all(x > y for x, y in zip(Vs_, Vs_[1:])) and all(x > p for x in Vs_) and Vs_[-1] - p < F(1, 10**6))
    check(f"p={p}: C_k strictly increasing, < c_crit = {cc}, c_crit - C_60 < 1e-6",
          all(x < y for x, y in zip(Cs, Cs[1:])) and all(x < cc for x in Cs) and cc - Cs[-1] < F(1, 10**6))

# optional: confirm the table cells of ../main.tex (if present) match
import os
for texname in ["../main.tex", "../csf_masks.tex"]:
    if os.path.exists(texname):
        tex = open(texname, encoding="utf8").read()
        cells = [x for p in TABLE for k in TABLE[p] for x in TABLE[p][k][2:] if k > 1]
        check(f"all 16 non-trivial table cells appear in {texname}", all(c in tex for c in cells))
        check(f"Example 14 values appear in {texname}",
              all(any(v in tex for v in (f"{n}/{d}", "{%s}{%s}" % (n, d)))
                  for n, d in [(3, 10), (79, 100), (39, 79), (58, 79), (177, 290), (4141, 5800)])
              and all(s in tex for s in ["0.4937", "0.7342", "0.6103", "0.7140", "\\ell^*=0", "a/(1-a)=3/4"]))
        break

# ------------------------------------------------------------------ Lemma 16 / Corollary 18 structure
print("== Lemma 16 increments, theta range, ties (Cor. 18)")
allok = True; tie_ok = True; notwo = True; trans_tie = True; theta_ok = True; E1 = True
grid = [F(n, 100) for n in range(51, 100)] + [F(2, 3), F(5, 8), F(4, 7)]
for p in grid:
    a = (1 - p) / p
    ls = lstar(a)
    E = {k: (Ck(a, k + 1) - Ck(a, k)) / (a**k / (1 - a**k) - a**(k + 1) / (1 - a**(k + 1))) for k in range(1, 40)}
    E1 &= E[1] == 1 - a
    ties = 0
    for k in range(1, 39):
        allok &= E[k] == Ek_formula(a, k)
        m = mk(a, k); tr = mk(a, k + 1) == m + 1
        tie = 2 * a**(m + 1) == 1 + a**k
        ties += tie
        if tie: trans_tie &= tr
        if tr: theta_ok &= 0 < theta(a, k) <= 1
        if k >= 2:
            trp = mk(a, k) == mk(a, k - 1) + 1
            notwo &= not (trp and tr)
            if trp: inc = (2 * a**m - 1 - a**k) + theta(a, k - 1) * (1 - a**(k - 1))
            elif tr: inc = (1 - theta(a, k)) * (1 - a**k)
            else: inc = 1 - a**k
            allok &= E[k] - E[k - 1] == inc
            allok &= (E[k] == E[k - 1]) == tie
    tie_ok &= ties <= ls
check("E_k from the C_k, V_k equals the formula of Thm 13 and the increments of Lemma 16 (52 biases, k<40)", allok)
check("E_1 = 1 - a", E1)
check("theta_k in (0,1] at every transition", theta_ok)
check("k-1 and k are never both transitions", notwo)
check("every tie is a transition; number of ties <= l*", trans_tie and tie_ok)

# ------------------------------------------------------------------ Corollary 21
print("== Corollary 21 (1/2 < p <= phi)")
c20 = True
for p in [F(51, 100), F(11, 20), F(3, 5), F(61, 100), F(618, 1000)]:
    a = (1 - p) / p
    c1 = (1 + a**2) / (1 + a + a**2); c2 = (1 + a**2 + 2 * a**3) / (1 + a + a**2 + a**3)
    c20 &= c1 == Ck(a, 3) and c2 == Ck(a, 4) and 1 - p < c1 < c2
    for c in [1 - p, (1 - p + c1) / 2, c1]:
        v = p + (1 - p)**2 / p * (1 - c)
        q = {-1: c - a * (1 - c), 0: 1 - c, 1: a * (1 - c)}
        lin = Vk_closed(p, 2) + (Vk_closed(p, 3) - Vk_closed(p, 2)) * (c - Ck(a, 2)) / (Ck(a, 3) - Ck(a, 2))
        c20 &= all(x >= 0 for x in q.values()) and sum(q.values()) == 1 and cost(q) == c and Vp(q, p) == v == lin
    for c in [c1, (c1 + c2) / 2, c2]:
        v = 1 - p + p**3 * (2 - c); t = p**2 * (2 - c)
        q = {-1: t, 0: a * t, 1: a**2 * t, 2: 1 - t * (1 + a + a**2)}
        lin = Vk_closed(p, 3) + (Vk_closed(p, 4) - Vk_closed(p, 3)) * (c - Ck(a, 3)) / (Ck(a, 4) - Ck(a, 3))
        c20 &= all(x >= 0 for x in q.values()) and sum(q.values()) == 1 and cost(q) == c and Vp(q, p) == v == lin
check("c_1 = C_3, c_2 = C_4, 1-p < c_1 < c_2; both affine pieces equal the interpolation; optimisers exact", c20)

# Dual-multiplier certificates of Corollary 21 (including tails n <= -2 and n >= 3)
c21_dual = True
for p in [F(51, 100), F(11, 20), F(3, 5), F(61, 100), F(618, 1000)]:
    a = (1 - p) / p
    # Piece 1: alpha_n = 1 (n<=0), alpha_1 = 1-a, alpha_n = 0 (n>=2)
    alpha1 = lambda n, a=a: F(1) if n <= 0 else (1 - a if n == 1 else F(0))
    for n in range(-10, 11):
        c21_dual &= 0 <= alpha1(n) <= 1
        kap = alpha1(n) + a * (1 - alpha1(n + 1))
        rhs = 1 + a**2 * (1 - abs(n))
        c21_dual &= (kap == rhs) if n in (-1, 0, 1) else (kap >= rhs)
    # Piece 2: alpha_n = 1 (n<=-1), alpha_0 = (2+a)p^2, alpha_1 = p^2, alpha_n = 0 (n>=2)
    alpha2 = lambda n, a=a, p=p: F(1) if n <= -1 else ((2 + a) * p**2 if n == 0 else (p**2 if n == 1 else F(0)))
    for n in range(-10, 11):
        c21_dual &= 0 <= alpha2(n) <= 1
        kap = alpha2(n) + a * (1 - alpha2(n + 1))
        rhs = a + p**2 * (2 - abs(n))
        c21_dual &= (kap == rhs) if n in (-1, 0, 1, 2) else (kap > rhs)
check("Corollary 21 dual multipliers alpha_n in [0,1] certify both affine global lower bounds", c21_dual)


# ------------------------------------------------------------------ Corollary 20 (tail of the frontier)
print("== Corollary 20 (tail: exact C_k, E_k for k >= k_0; k*(eps))")
def u_(a, k): return a**k / (1 - a**k)
def E_(p, k):
    a = (1 - p) / p
    return (Ck(a, k + 1) - Ck(a, k)) / (u_(a, k) - u_(a, k + 1))
def sigma_(p, k): return (2 * p - 1) / E_(p, k)
tail_ok = True; k0_ok = True; rows20 = []
for p in [F(52, 100), F(11, 20), F(3, 5), F(2, 3), F(7, 10), F(3, 4), F(4, 5), F(9, 10)]:
    a = (1 - p) / p; l = lstar(a); cc = ccrit(p)
    k0 = next(k for k in range(1, 10**4) if a**k < 2 * a**l - 1)
    kb = next(k for k in range(1, 10**4) if mk(a, k) == l)          # first k with m_k = l*
    k0_ok &= (k0 == kb) and all(mk(a, k) == l for k in range(k0, k0 + 25))
    for k in range(k0, k0 + 15):
        Cc = (cc - a**k * (k - l + a / (1 - a))) / (1 - a**k)
        cb, _, _ = brute(p, k)
        Ee = k - cc - l + a**(k + 1) / (1 - a)
        tail_ok &= Cc == Ck(a, k) == cb and E_(p, k) == Ee
        tail_ok &= 0 < E_(p, k) - (k - cc - l) <= a**(k + 1) / (1 - a)
    if k0 > 1:  # formula must fail just below k_0 (shows k_0 is sharp)
        k0_ok &= (cc - a**(k0 - 1) * (k0 - 1 - l + a / (1 - a))) / (1 - a**(k0 - 1)) != Ck(a, k0 - 1)
    rows20.append((p, l, k0))
check("k_0 = min{k: a^k < 2a^{l*}-1} equals the first k with m_k = l*, and is sharp: " +
      ", ".join(f"p={p}: l*={l}, k_0={k}" for p, l, k in rows20), k0_ok)
check("exact C_k (= brute force over shifts) and exact E_k = k - c_crit - l* + a^{k+1}/(1-a), k_0 <= k < k_0+15, 8 biases", tail_ok)
# k*(eps): brute force vs the corollary
kstar_ok = True; ex = {}
for p in [F(11, 20), F(3, 5), F(2, 3), F(7, 10), F(4, 5)]:
    a = (1 - p) / p; l = lstar(a); cc = ccrit(p)
    k0 = next(k for k in range(1, 10**4) if a**k < 2 * a**l - 1)
    sk0 = sigma_(p, k0)
    for eps in [F(1, 20), F(1, 10), F(1, 50), F(1, 100), F(1, 7), F(3, 100)]:
        if eps > sk0: continue
        kb = next(k for k in range(1, 10**4) if sigma_(p, k) < eps)
        x = (2 * p - 1) / eps + cc + l; fx = x.numerator // x.denominator
        kstar_ok &= x - a**(kb + 1) / (1 - a) < kb <= fx + 1
        if x - fx >= a**(fx + 1) / (1 - a): kstar_ok &= kb == fx + 1
        if eps == F(1, 20): ex[p] = kb
check("k*(eps) (brute force) satisfies the bounds of Corollary 20 for 6 values of eps and 5 biases", kstar_ok)
check(f"k*(0.05) = 7, 9, 13 for p = 0.6, 0.7, 0.8 (got {ex.get(F(3,5))}, {ex.get(F(7,10))}, {ex.get(F(4,5))})",
      (ex.get(F(3, 5)), ex.get(F(7, 10)), ex.get(F(4, 5))) == (7, 9, 13))

# ------------------------------------------------------------------ Fig. 3 caption, p = phi (exact, sympy)
print("== Fig. 3 caption, p = phi")
s5 = sp.sqrt(5)
pS = (s5 - 1) / 2; aS = sp.simplify((1 - pS) / pS)
z = lambda x: sp.simplify(sp.expand(x)) == 0
check("a = phi at p = phi", z(aS - pS))
check("tie: 2a = 1 + a^3 at p = phi (so k=3 is a tie with m_3 = 0)", z(2 * aS - 1 - aS**3))
def TS(L, k):
    g = (1 - aS) / (1 - aS**k); return {L + j: sp.nsimplify(sp.simplify(g * aS**j)) for j in range(k)}
def costS(q): return sp.simplify(sum(abs(n) * w for n, w in q.items()))
cs = {L: costS(TS(L, 3)) for L in range(-4, 2)}
cmin3 = min(cs.values(), key=lambda x: float(x))
mins3 = [L for L in cs if z(cs[L] - cmin3)]
cs2 = {L: costS(TS(L, 2)) for L in range(-3, 2)}
cmin2 = min(cs2.values(), key=lambda x: float(x))
mins2 = [L for L in cs2 if z(cs2[L] - cmin2)]
check(f"length-3 cost minimisers are L in {mins3} (expected [-1, 0])", sorted(mins3) == [-1, 0])
check(f"length-2 cost minimiser is L in {mins2} (expected [0])", mins2 == [0])
T02, T13 = TS(0, 2), TS(-1, 3)
mix = {n: sp.simplify((T02.get(n, 0) + T13.get(n, 0)) / 2) for n in (-1, 0, 1)}
shown = {-1: "0.250", 0: "0.464", 1: "0.286"}
check("masses of the equal mixture of T_{0,2} and T_{-1,3}: " + ", ".join(f"{float(mix[n]):.6f}" for n in (-1, 0, 1)),
      all("%.3f" % float(mix[n]) == shown[n] for n in mix) and z(sum(mix.values()) - 1))
def VpS(q):
    ys = set(q) | {y + 1 for y in q}
    tot = pS
    for y in ys:
        x = sp.expand((1 - pS) * q.get(y - 1, 0) - pS * q.get(y, 0))   # canonical r + s*sqrt(5)
        if x != 0 and float(x) > 0: tot += x
    return sp.expand(tot)
cmid = (cmin2 + cmin3) / 2
gmid = (pS + (2 * pS - 1) * (aS**2 / (1 - aS**2) + aS**3 / (1 - aS**3)) / 2)
check("the mixture has cost (C_2+C_3)/2 and V_p equal to the frontier value there (exact)",
      z(costS(mix) - cmid) and z(VpS(mix) - gmid))

# ------------------------------------------------------------------ c_crit asymptotics near p = 1/2
print("== c_crit near p = 1/2")
for p in [F(51, 100), F(501, 1000)]:
    r = float(ccrit(p)) / (math.log(2) / (2 * (2 * float(p) - 1)))
    print(f"   p={p}: c_crit={float(ccrit(p)):.4f}, ln2/(2(2p-1))={math.log(2)/(2*(2*float(p)-1)):.4f}, ratio {r:.4f}")

n_fail = sum(1 for _, ok in RESULTS if not ok)
print(f"\nverify_numbers.py: {len(RESULTS) - n_fail}/{len(RESULTS)} checks passed -> " + ("ALL PASS" if n_fail == 0 else "FAIL"))
raise SystemExit(1 if n_fail else 0)
