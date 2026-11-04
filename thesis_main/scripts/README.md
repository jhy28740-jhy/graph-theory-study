# Numerical scripts for Chapter 7

All scripts are deterministic and depend only on `numpy` (plus the standard
library). They reproduce every number quoted in Chapter 7 and in the
"Additions made in this thesis" part of Appendix A.

| script | produces | thesis location |
|---|---|---|
| `verify_constants_exact.py` | 60-digit values of `c`, `c_w`, `(3+sqrt5)/4`, `3/sqrt5`, the equal-blow-up coefficient, `b*`, `p*`, and the check `Sigma_3(quotient at b*) == c_w` | Section 7.1, Appendix A |
| `exhaustive_n3_n6.py` | all 33864 labelled graphs for `3 <= n <= 6`, extrema and isomorphism classes | Table 7.1 |
| `exhaustive_n7.py` | all 2097152 labelled graphs on 7 vertices, extremum and structure of the 252 maximisers | Table 7.1 |
| `blowup_families.py` | equal and optimally weighted blow-ups of `K1 v C5` | Tables 7.2, 7.3 |
| `weighted_random.py` | randomised test of the weighted form, plus the negative control | Table 7.4 |

Raw outputs are written to `../data/*.json`.

These computations are **numerical evidence only**. They are not used as a
step in any proof; see Section 7.5 of the thesis.
