# Optimal Additive Integer Masks for One-Guess Vulnerability under an Absolute-Moment Budget

Joonas Pääkkönen, Dalarna University

Manuscript, figure and verification code for the paper. Version 2 (release v1.1.0) determines the complete cost–vulnerability frontier for a biased-bit secret (Theorem 5.8). All curves in the frontier figure are proved.

arXiv: to be added once announced.

## Contents
- `manuscript/`: `mask_arxiv_manuscript.tex`, `mask_arxiv_manuscript.pdf`, `frontier.pdf`
- `examples/`: `toy-example.tex`, `toy-example.pdf`
- `code/`: v1 scripts (`lp_upper.py`, `check_certs.py`, `simulate.py`); the dual certificates checked by `check_certs.py` remain valid in v2
- `code/v2/`: verification scripts for v2 (see Remark 5.19)
  - `frontier_exact.py`: exact rational points Q_k = (C_k, V_k) and the frontier
  - `verify_theory.py`: exact-arithmetic checks of the algebraic ingredients of Theorem 5.8; writes `data/verify_theory.json`
  - `core.py`: LP and certificate routines (integer and modular)
  - `run_numerics.py`: floating-point LP cross-checks on a (p, c) grid; writes `data/numerics_grid.json`
  - `modular_n_needed.py`: smallest modulus whose modular LP matches the frontier (numerical evidence only); writes `data/modular_n_needed.json`
  - `make_frontier_v2.py`: regenerates `frontier.pdf`
- `LICENSE` (MIT, code) and `LICENSE-manuscript` (CC BY 4.0, manuscript and example)

## Reproduce
Requires Python 3 with numpy, scipy and matplotlib.

    cd code/v2
    python3 verify_theory.py
    python3 make_frontier_v2.py

## Compile
    cd manuscript && pdflatex mask_arxiv_manuscript.tex

## License
Code: MIT. Manuscript text: CC BY 4.0.
