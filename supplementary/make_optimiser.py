"""Regenerate optimiser.pdf (Fig. 3 of the CSF version): the optimal mask at the midpoint of the piece
[C_2, C_3] for p = phi = (sqrt5 - 1)/2, i.e. the equal mixture of T_{0,2} and T_{-1,3}.
At p = phi one has a = phi and 2a = 1 + a^3, so k = 3 is a tie and the shift -1 is as cheap as -m_3 = 0
(Lemma 15(2)).  The masses are computed exactly in Q(sqrt 5) (sympy) and checked in verify_numbers.py.
Run from this folder: python3 make_optimiser.py   (writes ../optimiser.pdf)
"""
import sympy as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

phi = (sp.sqrt(5) - 1) / 2
a = sp.simplify((1 - phi) / phi)            # = phi
def T(L, k):
    g = (1 - a) / (1 - a**k)
    return {L + j: sp.simplify(g * a**j) for j in range(k)}
T02, T13 = T(0, 2), T(-1, 3)
q = {n: sp.simplify((T02.get(n, 0) + T13.get(n, 0)) / 2) for n in (-1, 0, 1)}

if __name__ == "__main__":
    for n in q: print(n, q[n], float(q[n]))
    fig, ax = plt.subplots(figsize=(4.0, 2.6))
    ns = sorted(q)
    ax.bar(ns, [float(q[n]) for n in ns], color="steelblue", edgecolor="black", width=0.8)
    for n in ns:
        ax.text(n, float(q[n]) + 0.01, "%.3f" % float(q[n]), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(ns); ax.set_xlabel(r"$n$"); ax.set_ylabel(r"$q_c(n)$")
    ax.set_ylim(0, 0.55)
    ax.set_title(r"Optimal mask on the piece $[C_2,C_3]$ ($p=\varphi$)", fontsize=10)
    fig.tight_layout()
    meta = {"Creator": None, "Producer": None, "CreationDate": None}
    fig.savefig("../optimiser.pdf", metadata=meta)
