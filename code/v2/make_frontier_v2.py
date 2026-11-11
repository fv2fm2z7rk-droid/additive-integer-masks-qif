"""Regenerate the biased-bit frontier figure for v2: the full proven frontier
(Theorem 'complete frontier') is drawn solid, with the vertices Q_k as dots.
Exact rational arithmetic via frontier_exact.py (same folder)."""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction as F
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from frontier_exact import envelope_points, H, ccrit

PS = [(F(11, 20), "p=0.55"), (F(3, 5), "p=0.60"), (F(2, 3), "p=2/3"), (F(3, 4), "p=0.75"), (F(9, 10), "p=0.90")]
CMAX = 3.6
fig, ax = plt.subplots(figsize=(6.4, 3.6))
for (pF, lab), col in zip(PS, ["C0", "C1", "C2", "C3", "C4"]):
    hull = envelope_points(pF, 300)
    cs = np.linspace(0, CMAX, 1500)
    ax.plot(cs, [float(H(pF, F(c).limit_denominator(10**9), hull=hull)) for c in cs], "-", color=col, lw=1.4, label=lab)
    pts = [(float(x), float(y)) for x, y, k in hull if k != "inf" and float(y) - float(pF) > 1e-4]
    ax.plot([x for x, _ in pts], [y for _, y in pts], "o", color=col, ms=2.5)
    ax.plot([float(ccrit(pF))], [float(pF)], "|", color=col, ms=9, mew=1.5)
ax.set_xlim(0, CMAX); ax.set_ylim(0.5, 1.03)
ax.set_xlabel(r"$c=\mathbb{E}|Y|$"); ax.set_ylabel(r"$V_p^*(c)$")
ax.set_title("Cost--vulnerability frontier (all curves proved)\ndots: vertices $Q_k$; bars: $c_{\\mathrm{crit}}(p)$", fontsize=9)
ax.legend(fontsize=8)
fig.tight_layout(); fig.savefig("frontier.pdf"); fig.savefig("frontier_preview.png", dpi=150)
