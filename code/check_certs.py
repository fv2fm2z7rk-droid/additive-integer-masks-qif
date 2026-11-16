#!/usr/bin/env python3
"""Exact-arithmetic checks for the two intermediate certificates and examples.
Does not solve an LP. Does not claim to prove the certificates for all parameters.
"""
from fractions import Fraction

def check_p35():
    p = Fraction(3, 5)
    a = (1 - p) / p
    c1 = (1 + a**2) / (1 + a + a**2)
    c2 = (1 + a**2 + 2*a**3) / (1 + a + a**2 + a**3)
    assert c1 == Fraction(13, 19)
    assert c2 == Fraction(11, 13)
    # V at c=1/2
    assert Fraction(13, 15) - Fraction(4, 15)*Fraction(1, 2) == Fraction(11, 15)
    # V at c=3/4
    assert Fraction(2, 5) + Fraction(27, 125)*(2 - Fraction(3, 4)) == Fraction(67, 100)
    # V at c=4/5
    assert Fraction(2, 5) + Fraction(27, 125)*(2 - Fraction(4, 5)) == Fraction(412, 625)
    print("p=3/5 exact identities: OK")

def check_counterexample():
    p = Fraction(3, 4)
    c = Fraction(3, 10)
    formula = (p**2 + (1 - c)*(1 - p)**2) / p
    # mask (51/70, 17/70, 2/70)
    q0, q1, q2 = Fraction(51, 70), Fraction(17, 70), Fraction(2, 70)
    V = p*q0 + max(p*q1, (1-p)*q0) + max(p*q2, (1-p)*q1) + (1-p)*q2
    assert V < formula
    print(f"counterexample: formula={float(formula):.6f}, mask V={float(V):.6f}")

def check_dual_conditions():
    # first certificate requires (1-p)^2 >= p*(2p-1) i.e. p <= phi
    # second requires p(1+p) <= 1 i.e. p <= phi
    print("dual feasibility for the two certificates holds precisely for p <= (sqrt(5)-1)/2")
    print("(analytic fact from the paper; not re-proved by finite testing)")

if __name__ == "__main__":
    check_p35()
    check_counterexample()
    check_dual_conditions()
    print("checks completed (exact identities and counterexample only)")
