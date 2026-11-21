"""Regenerate ccrit.pdf (Fig. 2 of the CSF version): critical cost c_crit(p) versus bias p.
c_crit(p) = min_{l>=0} ( l + a(2a^l - 1)/(1-a) ),  a = (1-p)/p   (Theorem 9).
"""
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def ccrit(p):
    a = (1 - p) / p
    lstar = math.ceil(math.log(2) / math.log(1 / a)) - 1   # integer median (Theorem 9)
    best = min(l + a * (2 * a**l - 1) / (1 - a) for l in range(max(0, lstar - 2), lstar + 3))
    return best

if __name__ == "__main__":
    P = np.linspace(0.51, 0.99, 600)
    C = [ccrit(p) for p in P]
    fig, ax = plt.subplots(figsize=(4.0, 2.8))
    ax.plot(P, C, color="blue", lw=2)
    ax.set_xlabel(r"$p$")
    ax.set_ylabel(r"$c_{\mathrm{crit}}(p)$")
    ax.set_title("Critical cost versus bias")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig("../ccrit.pdf", bbox_inches="tight", pad_inches=0.05,
                metadata={"Creator": None, "Producer": None, "CreationDate": None})
    for p in [0.51, 0.501, 0.5001, 0.6, 2/3, 0.75, 0.9]:
        print(p, ccrit(p), math.log(2) / (2 * (2 * p - 1)))
