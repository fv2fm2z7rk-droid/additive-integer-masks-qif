# Verification Scripts: Optimal Additive Integer Masks for One-Guess Vulnerability

Exact-arithmetic verification scripts, weak-duality LP checks, and figure-generation code for:

> **Joonas Pääkkönen**, *Optimal Additive Integer Masks for One-Guess Vulnerability under an Absolute-Moment Budget*, arXiv:2611.16066 (v3, Nov. 2026).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24418378.svg)](https://doi.org/10.5281/zenodo.24418378)

## Requirements

- Python 3
- `numpy`, `scipy`, `sympy`, `matplotlib`

Install dependencies with:

```bash
pip install -r requirements.txt
```

## Quick Start

All scripts use relative imports and should be run from the `supplementary/` directory:

```bash
cd supplementary
python3 verify_numbers.py
python3 dual_grid_check.py
python3 levelset_check.py
python3 check_review_fixes.py
python3 check_peeling_gap.py
```

To regenerate Figs. 1–3 of the paper (`frontier.pdf`, `ccrit.pdf`, `optimiser.pdf` in the repository root):

```bash
cd supplementary
python3 make_frontier.py
python3 make_ccrit.py
python3 make_optimiser.py
```

## Contents

| File | Purpose |
|---|---|
| `supplementary/verify_numbers.py` | Recomputes every hand-typed number in the paper (Example 14, Appendix B table, Fig. 3 exact $\mathbb{Q}(\sqrt{5})$ masses, Corollary 20 tail bounds and $k^*(\varepsilon)$, Corollary 21 endpoints, optimisers, and dual certificates $\alpha_n \in [0,1]$) in exact arithmetic and against brute force over shifts; also checks Lemma 16 increments and the Corollary 18 tie bound. |
| `supplementary/dual_grid_check.py` | Compares rigorous rationalised weak-duality LP lower bounds (all $y \in \mathbb{Z}$, tails included) and primal LP upper bounds against the Theorem 13 closed-form frontier on a 64-point $(p,c)$ grid. |
| `supplementary/core.py` | Primal/dual LP and exact rational certificate routines used by `dual_grid_check.py`. |
| `supplementary/levelset_check.py` | Exact-arithmetic checks of the level-set layer-cake decomposition (Lemma A) on random finitely supported masks and long two-sided-tail masks. |
| `supplementary/check_review_fixes.py` | Exact checks of Lemma 15(1)–(2) and shared helper definitions (`T`, `Vp`, `cost`, `mk`, `A`, `Ck`). |
| `supplementary/check_peeling_gap.py` | Explicit counterexample showing greedy peeling is not additive for $V_p$ (motivating Lemma A). |
| `supplementary/make_frontier.py`, `make_ccrit.py`, `make_optimiser.py` | Regenerate Figs. 1–3 with stripped PDF metadata. |

## Citation

```bibtex
@misc{Paakkonen2026masks,
  author       = {P{\"a}{\"a}kk{\"o}nen, Joonas},
  title        = {Optimal Additive Integer Masks for One-Guess Vulnerability under an Absolute-Moment Budget},
  year         = {2026},
  eprint       = {2611.16066},
  archivePrefix = {arXiv}
}
```

## License

MIT License. See [LICENSE](LICENSE).
