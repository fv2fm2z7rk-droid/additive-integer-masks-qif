# Supplementary scripts

All scripts use only relative imports and are run from this folder (`python3 <script>.py`).
Requirements: Python 3, sympy, numpy, scipy, matplotlib. Figure scripts write into the parent folder.

| Script | Purpose |
|---|---|
| `verify_numbers.py` | Recomputes every hand-typed number in the paper (Example 14, appendix table, Fig. 3 masses, Corollary 21) in exact arithmetic and against brute force over shifts; also checks the dual certificates of Corollary 21, Lemma 16 increments, the tie bound of Corollary 18, and the tail formulas and k*(eps) of Corollary 20. Prints PASS/FAIL. |
| `check_review_fixes.py` | Exact checks of Lemma 15(1)–(2); definitions `T`, `Vp`, `cost`, `mk`, `A`, `Ck` shared by the other scripts. Part (c) shows that a greedy peeling is not additive for `V_p`. |
| `check_peeling_gap.py` | Counterexample to additivity of the greedy peeling (the reason for Lemma A). |
| `levelset_check.py` | Exact checks of the level-set layer cake of Lemma A on random masks, including long two-sided tails. |
| `dual_grid_check.py` | LP lower (weak-duality certificate) and upper bounds compared with the closed-form frontier of Theorem 13 on a 64-point grid; uses `core.py`. |
| `core.py` | LP and certificate routines. |
| `make_frontier.py`, `make_ccrit.py`, `make_optimiser.py` | Regenerate Figs. 1–3 (`../frontier.pdf`, `../ccrit.pdf`, `../optimiser.pdf`). |
