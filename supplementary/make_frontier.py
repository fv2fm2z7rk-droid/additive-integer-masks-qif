"""Regenerate frontier.pdf (Fig. 1 of the CSF version): the biased-bit frontier V_p^*(c) of Theorem 13,
i.e. the linear interpolation of the points Q_k = (C_k, V_k) of Definition 12, constant p for c >= c_crit(p).
Dots: vertices Q_k (drawn while V_k - p > 1e-4); vertical ticks: (c_crit(p), p).
Run from this folder: python3 make_frontier.py   (writes ../frontier.pdf)
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fractions import Fraction as F
from check_review_fixes import Ck

def Vk(p, k):
    a = (1 - p) / p
    return p + (2 * p - 1) * a**k / (1 - a**k)

def ccrit(p):
    a = (1 - p) / p
    l = 0
    while a**(l + 1) > 0.5: l += 1          # l* = #{j >= 1 : a^j > 1/2}  (Theorem 9)
    return l + a * (2 * a**l - 1) / (1 - a)

PS = [(0.55, "$p=0.55$"), (0.60, "$p=0.60$"), (2 / 3, "$p=2/3$"), (0.75, "$p=0.75$"), (0.90, "$p=0.90$")]
CMAX = 3.6
if __name__ == "__main__":
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    for (p, lab), col in zip(PS, ["C0", "C1", "C2", "C3", "C4"]):
        a = (1 - p) / p
        pts = [(Ck(a, k), Vk(p, k)) for k in range(1, 400)]
        cc = ccrit(p)
        xs = [x for x, _ in pts] + [cc, CMAX]
        ys = [y for _, y in pts] + [p, p]
        ax.plot(xs, ys, "-", color=col, lw=1.4, label=lab)
        dots = [(x, y) for x, y in pts if y - p > 1e-4]
        ax.plot([x for x, _ in dots], [y for _, y in dots], "o", color=col, ms=2.5)
        ax.plot([cc], [p], "|", color=col, ms=9, mew=1.5)
        print(f"p={p:.4f}: c_crit={cc:.6f}, C_399={pts[-1][0]:.6f}, V_399-p={pts[-1][1]-p:.2e}")
    ax.set_xlim(0, CMAX); ax.set_ylim(0.5, 1.02)
    ax.set_xlabel(r"budget $c=\mathbb{E}|N|$"); ax.set_ylabel(r"$V_p^*(c)$")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0, frameon=False)
    fig.tight_layout()
    meta = {"Creator": None, "Producer": None, "CreationDate": None}
    fig.savefig("../frontier.pdf", metadata=meta)
