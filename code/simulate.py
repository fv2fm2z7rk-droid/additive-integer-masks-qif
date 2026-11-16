import numpy as np

np.random.seed(0)
n = 100000
p = 0.6

def simulate(q, support):
    Y = np.random.choice(support, size=n, p=q)
    X = np.random.choice([0, 1], size=n, p=[p, 1-p])
    S = X + Y
    success = 0
    for s in np.unique(S):
        mask = (S == s)
        p0 = np.mean(X[mask] == 0)
        guess = 0 if p0 >= 0.5 else 1
        success += np.sum(X[mask] == guess)
    return success / n

q_bad = [0.5, 0.5]
supp_bad = [0, 1]
v_bad = simulate(q_bad, supp_bad)

q_good = [1/6, 1/2, 1/3]
supp_good = [-1, 0, 1]
v_good = simulate(q_good, supp_good)

print(f"Bad mask:  estimated V = {v_bad:.4f} (exact 0.80)")
print(f"Good mask: estimated V = {v_good:.4f} (exact 11/15 ≈ 0.7333)")
