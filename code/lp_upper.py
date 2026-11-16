#!/usr/bin/env python3
"""Finite-support LP upper bounds for the biased-bit problem.

Correct epigraph constraints:
  p*q_s - z_s <= 0
  (1-p)*q_(s-1) - z_s <= 0

Requires: scipy, numpy.
A truncated optimum is an upper bound only; it does not certify the unrestricted minimum.
"""
import numpy as np
from scipy.optimize import linprog

def solve_lp(p, c, M=40, tol=1e-8):
    nq = 2*M + 1
    nz = 2*M + 2
    n = nq + nz
    obj = np.zeros(n)
    obj[nq:] = 1.0

    A_eq = [np.zeros(n)]
    A_eq[0][:nq] = 1.0
    b_eq = [1.0]

    A_ub = []
    b_ub = []
    row = np.zeros(n)
    for i, y in enumerate(range(-M, M+1)):
        row[i] = abs(y)
    A_ub.append(row)
    b_ub.append(c)

    for s_idx, s in enumerate(range(-M, M+2)):
        row = np.zeros(n)
        if -M <= s <= M:
            row[s + M] = p
        row[nq + s_idx] = -1.0
        A_ub.append(row)
        b_ub.append(0.0)
        row = np.zeros(n)
        if -M <= s-1 <= M:
            row[(s-1) + M] = (1-p)
        row[nq + s_idx] = -1.0
        A_ub.append(row)
        b_ub.append(0.0)

    bounds = [(0, None)] * n
    res = linprog(obj, A_ub=np.array(A_ub), b_ub=np.array(b_ub),
                  A_eq=np.array(A_eq), b_eq=np.array(b_eq),
                  bounds=bounds, method='highs')

    if not res.success:
        raise RuntimeError(f"LP failed: {res.message}")
    if not np.all(np.isfinite(res.x)):
        raise RuntimeError("non-finite solution")

    q = res.x[:nq]
    V_lp = res.fun

    if abs(q.sum() - 1.0) > tol:
        raise RuntimeError(f"normalization failed: sum q = {q.sum()}")
    if np.any(q < -tol):
        raise RuntimeError("negative mass")
    cost = sum(abs(y)*q[i] for i, y in enumerate(range(-M, M+1)))
    if cost > c + tol:
        raise RuntimeError(f"budget violated: {cost} > {c}")

    V2 = 0.0
    for s in range(-M, M+2):
        qs = q[s+M] if -M <= s <= M else 0.0
        qsm = q[s-1+M] if -M <= s-1 <= M else 0.0
        V2 += max(p*qs, (1-p)*qsm)
    if abs(V_lp - V2) > 1e-6:
        raise RuntimeError(f"objective mismatch: LP={V_lp}, recomputed={V2}")
    if not (p - tol <= V2 <= 1 + tol):
        raise RuntimeError(f"vulnerability out of range: {V2}")

    return V_lp, V2, q

def check(p, c, reference, tol, M=40, label=""):
    V, V2, q = solve_lp(p, c, M=M)
    if abs(V - reference) > tol:
        raise RuntimeError(
            f"regression failed ({label}): p={p}, c={c}, M={M}, "
            f"actual={V:.8f}, reference={reference:.8f}, tol={tol}"
        )
    print(f"OK  p={p:.4f}, c={c:.2f}, V={V:.6f}")

if __name__ == "__main__":
    # numerical comparison values (rounded to 6 decimals; tolerance accommodates rounding)
    check(3/5, 1.0, 0.631429, 5e-6, label="p=3/5 c=1")
    check(3/5, 1.5, 0.603591, 5e-6, label="p=3/5 c=1.5")
    check(2/3, 0.5, 0.733333, 5e-6, label="p=2/3 c=0.5")
    check(2/3, 0.75, 0.687075, 5e-6, label="p=2/3 c=0.75")
    # analytic exact value
    check(3/5, 0.75, 0.67, 1e-8, label="p=3/5 c=0.75 (exact 67/100)")
    # additional checks
    check(0.6, 0.0, 1.0, 1e-8, label="c=0")
    check(0.6, 0.2, 1 - 0.6*0.2, 1e-6, label="initial segment")
    print("all regression tests passed")
