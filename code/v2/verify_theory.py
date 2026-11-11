"""Exact (rational-arithmetic) checks of every algebraic ingredient of the
gap-closure theorem.  Writes data/verify_theory.json."""
import json, random
from fractions import Fraction as F
from frontier_exact import Qk, Ck_bruteforce, ccrit, envelope_points, H, m_of
from core import exact_V

random.seed(1)
P_GRID = [F(j, 200) for j in range(101, 200)] + [F(618034, 1000000), F(618033, 1000000), F(2, 3)]
out = {}

def T(a, L, k):
    t = (1 - a) / (1 - a ** k)
    return {L + j: t * a ** j for j in range(k)}

def Ck_closed(p, k):
    a = (1 - p) / p; m = m_of(a, k)
    u = a ** k / (1 - a ** k); A = a * (1 - a ** m) / (1 - a)
    return a / (1 - a) + (m - 2 * A) + u * (2 * (m - A) - k)

# 1) C_k closed form, V_k formula, minimiser m_k
bad = 0; n = 0
for p in P_GRID:
    a = (1 - p) / p
    for k in range(1, 26):
        C, V, m = Qk(p, k)
        n += 1
        if C != Ck_bruteforce(p, k) or C != Ck_closed(p, k): bad += 1
        q = T(a, -m, k)
        if exact_V(p, q) != V or V != p + (2 * p - 1) * a ** k / (1 - a ** k): bad += 1
        if sum(abs(y) * v for y, v in q.items()) != C: bad += 1
out['closed_forms_checked'] = n; out['closed_form_failures'] = bad

# 2) E_k = (C_{k+1}-C_k)/(u_k-u_{k+1}) nondecreasing, E_1 = 1-a, no consecutive transitions,
#    equality E_k = E_{k-1} only at ties 2a^{m+1} = 1+a^k
viol = 0; consec = 0; ties = []; n = 0
for p in P_GRID:
    a = (1 - p) / p; K = 60 if p > F(56, 100) else 40
    u = lambda k: a ** k / (1 - a ** k)
    Cs = {k: Qk(p, k)[0] for k in range(1, K + 2)}
    ms = {k: m_of(a, k) for k in range(1, K + 2)}
    E = {k: (Cs[k + 1] - Cs[k]) / (u(k) - u(k + 1)) for k in range(1, K + 1)}
    if E[1] != 1 - a: viol += 1
    for k in range(2, K + 1):
        n += 1
        if E[k] < E[k - 1]: viol += 1
        if E[k] == E[k - 1]: ties.append((str(p), k, 2 * a ** (ms[k] + 1) == 1 + a ** k))
    for k in range(1, K):
        if ms[k + 1] - ms[k] > 1: viol += 1
        if ms[k + 1] > ms[k] and ms[k + 2] > ms[k + 1]: consec += 1
out['E_monotone_checks'] = n; out['E_monotone_violations'] = viol
out['consecutive_transitions'] = consec; out['equalities_found'] = ties

# 3) layer-cake decomposition identity on random rational masks
def decompose(p, q):
    a = (1 - p) / p
    w = {s: v / a ** s for s, v in q.items() if v > 0}
    levels = sorted(set(w.values()))
    pieces = []   # (weight dt, L, k)
    prev = F(0); S = sorted(w)
    for lev in levels:
        dt = lev - prev; prev = lev
        sup = [s for s in S if w[s] >= lev]
        # components of {w >= lev} = {w > t} for t in [prev_old, lev)
        comps = []; start = sup[0]; last = sup[0]
        for s in sup[1:]:
            if s == last + 1: last = s
            else: comps.append((start, last)); start = last = s
        comps.append((start, last))
        for L, R in comps: pieces.append((dt, L, R - L + 1))
    return pieces
fails = 0; hull_viol = 0; trials = 0
for p in [F(11, 20), F(3, 5), F(13, 20), F(3, 4), F(9, 10)]:
    a = (1 - p) / p; hull = envelope_points(p, 80)
    for _ in range(150):
        lo = random.randint(-6, 2); hi = lo + random.randint(0, 9)
        q = {s: F(random.randint(0, 20)) for s in range(lo, hi + 1)}
        tot = sum(q.values())
        if tot == 0: continue
        q = {s: v / tot for s, v in q.items() if v}
        pcs = decompose(p, q)
        mass = sum(dt * a ** L * (1 - a ** k) / (1 - a) for dt, L, k in pcs)
        Vmix = sum(dt * a ** L * (1 - a ** k) / (1 - a) * (p + (2 * p - 1) * a ** k / (1 - a ** k)) for dt, L, k in pcs)
        cost = sum(dt * sum(abs(L + j) * a ** (L + j) for j in range(k)) for dt, L, k in pcs)
        recon = {}
        for dt, L, k in pcs:
            for j in range(k): recon[L + j] = recon.get(L + j, 0) + dt * a ** (L + j)
        trials += 1
        if mass != 1 or Vmix != exact_V(p, q) or cost != sum(abs(s) * v for s, v in q.items()) or recon != q:
            fails += 1
        if exact_V(p, q) < H(p, sum(abs(s) * v for s, v in q.items()), hull=hull): hull_viol += 1
out['decomposition_trials'] = trials; out['decomposition_failures'] = fails; out['masks_below_envelope'] = hull_viol

# 4) agreement with the integer paper (I1)-(I3) and its numerical examples
chk = {}
for p in [F(52, 100), F(55, 100), F(3, 5), F(61, 100)]:
    a = (1 - p) / p; hull = envelope_points(p, 80)
    c1 = (1 + a * a) / (1 + a + a * a); c2 = (1 + a * a + 2 * a ** 3) / (1 + a + a * a + a ** 3)
    ok = True
    for i in range(11):
        c = (1 - p) + (c1 - (1 - p)) * F(i, 10)
        ok &= H(p, c, hull=hull) == p + (1 - p) ** 2 / p * (1 - c)
        c = c1 + (c2 - c1) * F(i, 10)
        ok &= H(p, c, hull=hull) == 1 - p + p ** 3 * (2 - c)
        c = (1 - p) * F(i, 10)
        ok &= H(p, c, hull=hull) == 1 - p * c
    ok &= Qk(p, 3)[0] == c1 and Qk(p, 4)[0] == c2
    chk[str(p)] = ok
for p in [F(3, 4), F(9, 10), F(2, 3)]:
    hull = envelope_points(p, 80)
    chk[str(p)] = all(H(p, (1 - p) * F(i, 10), hull=hull) == 1 - p * (1 - p) * F(i, 10) for i in range(11))
out['I1_I3_agreement'] = chk
ex = {}
for p, c in [(F(3, 5), F(1, 2)), (F(3, 5), F(3, 4)), (F(3, 5), F(4, 5)), (F(3, 5), F(1)), (F(3, 5), F(3, 2)),
             (F(2, 3), F(1, 2)), (F(2, 3), F(3, 4)), (F(3, 4), F(3, 10)), (F(3, 4), F(1, 4)), (F(11, 20), F(3, 2))]:
    hull = envelope_points(p, 120)
    v = H(p, c, hull=hull)
    k = max(h[2] for h in hull if h[2] != 'inf' and h[0] <= c)
    ex[f"p={p},c={c}"] = dict(value=str(v), float=float(v), cell_k=k)
out['examples'] = ex
# p = 3/4, c = 3/10: paper's feasible mask value 223/280 vs. optimum
out['p34_c310_mask_value'] = str(exact_V(F(3, 4), {0: F(51, 70), 1: F(17, 70), 2: F(2, 70)}))
# c_crit consistency: limit of C_k
out['ccrit_check'] = {str(p): float(ccrit(p) - Qk(p, 60)[0]) for p in [F(3, 5), F(3, 4), F(11, 20)]}
json.dump(out, open('data/verify_theory.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
