"""Exact (rational) truncated-geometric points Q_k=(C_k,V_k) and the envelope H_p."""
from fractions import Fraction as F
import math

def m_of(a, k):
    return sum(1 for j in range(1, k) if 2 * a ** j > 1 + a ** k)

def Qk(p, k):
    """Return (C_k, V_k, m_k) for T_{-m_k,k}; works for Fractions or floats."""
    a = (1 - p) / p
    m = m_of(a, k)
    t = (1 - a) / (1 - a ** k)
    C = t * sum(abs(j - m) * a ** j for j in range(k))
    V = p + p * (1 - a) * a ** k / (1 - a ** k)
    return C, V, m

def Ck_bruteforce(p, k):
    a = (1 - p) / p; t = (1 - a) / (1 - a ** k)
    return min(t * sum(abs(L + j) * a ** j for j in range(k)) for L in range(-k, 1))

def ccrit(p):
    a = (1 - p) / p
    best = None; l = 0
    while True:
        v = l + a * (2 * a ** l - 1) / (1 - a)
        if best is None or v < best: best = v
        elif 2 * a ** (l + 1) < 1: break
        l += 1
    return best

def lower_hull(pts):
    """pts: list of (x,y,label); returns lower-left monotone hull (exact orientation)."""
    P = sorted(pts, key=lambda z: (z[0], z[1])); H = []
    for z in P:
        while len(H) >= 2:
            (x1, y1, _), (x2, y2, _) = H[-2], H[-1]
            if (x2 - x1) * (z[1] - y1) - (y2 - y1) * (z[0] - x1) <= 0: H.pop()
            else: break
        H.append(z)
    # keep only nonincreasing part (y decreasing in x)
    out = [H[0]]
    for z in H[1:]:
        if z[1] < out[-1][1]: out.append(z)
    return out

def envelope_points(p, K):
    pts = [(*Qk(p, k)[:2], k) for k in range(1, K + 1)]
    pts.append((ccrit(p), p, 'inf'))
    return lower_hull(pts)

def H(p, c, K=None, hull=None):
    if hull is None: hull = envelope_points(p, K)
    if c >= hull[-1][0]: return hull[-1][1]
    for (x1, y1, k1), (x2, y2, k2) in zip(hull, hull[1:]):
        if x1 <= c <= x2:
            return y1 + (y2 - y1) * (c - x1) / (x2 - x1)
    raise ValueError
