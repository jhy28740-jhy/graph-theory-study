# Numbering cross-table: preprint -> thesis

This file records the provenance of every numbered result in the thesis
`thesis_main/main.tex` relative to the preprint it expands,
`paper/top3_joint_moment_bound.tex` (H. Ji, *A joint moment bound for the three
largest graph eigenvalues*, 3 November 2026, 22 pages).

How to read it. The preprint column uses the preprint's own theorem/section
numbers and page numbers; the thesis column uses the thesis numbers (theorems
are numbered per chapter) and the page numbers of the compiled `main.pdf`
(60 pages). The `label` column gives the LaTeX label, which is the canonical
way to cite a result: no theorem number is ever hard-coded in the thesis prose.

Structure of the table.
- **Part 1** lists the 24 results that are carried over from the preprint. These
  blocks are extracted verbatim by line range into `thesis_main/ported/*.tex`
  (each ported file records its source line range in a header comment); the only
  edits are cross-reference adaptations such as `\S6` -> `Section~\ref{sec-quartic}`
  and `\S\ref{sec-further}` -> `Chapter~\ref{sec-further}`, one display-math line
  split in `p_moments_reg.tex`, and the restoration of a lead-in sentence in
  `p_intro_results.tex`.
- **Part 2** maps preprint sections to thesis chapters.
- **Part 3** lists labels that exist only in the thesis, i.e. genuinely new
  material: the literature review of Chapter 1, the preliminaries of Chapter 2
  (Ky Fan, graphons, blow-ups, harmonic kernels, and the projection/isotropic-law
  bridge Lemma 2.3), the weighted `[0,1]` scope results of Section 6.4, and the
  whole of Chapter 7 (numerical verification protocol, Tables 7.1--7.4, and
  Proposition 7.1 on equal blow-ups).

Reproduction. The thesis page numbers and numbers below were read from a fresh
three-pass `pdflatex` build of `thesis_main/main.tex`; the preprint numbers were
read from a fresh build of `paper/top3_joint_moment_bound.tex`. The auxiliary
scripts used for the Chapter 7 numerics are in `thesis_main/scripts/`, with
their JSON outputs committed in `thesis_main/data/`.

| preprint number | preprint page | thesis number | thesis page | label | kind |
|---|---|---|---|---|---|
| 1.1 | 2 | 1.1 | 3 | `main` | result |
| 1.2 | 2 | 1.2 | 4 | `structure-main` | result |
| 2.1 | 3 | 3.1 | 13 | `ranktwo` | result |
| 2.2 | 4 | 3.2 | 14 | `attainment` | result |
| 2.3 | 5 | 3.3 | 14 | `rem-newton` | result |
| 2.4 | 5 | 3.4 | 15 | `ellipse-structure` | result |
| 2.5 | 5 | 3.5 | 15 | `row-insertion` | result |
| 2.6 | 6 | 3.6 | 16 | `contact-structure` | result |
| 3.1 | 7 | 4.1 | 19 | `negative` | result |
| 3.2 | 7 | 4.2 | 19 | `propagate` | result |
| 3.3 | 8 | 4.3 | 20 | `descent` | result |
| 3.4 | 8 | 4.4 | 20 | `finite` | result |
| 4.2 | 9 | 4.6 | 22 | `binary` | result |
| 4.3 | 10 | 4.7 | 23 | `six` | result |
| 4.4 | 10 | 4.8 | 23 | `twentyfour` | result |
| 6.1 | 13 | 5.1 | 27 | `scalar` | result |
| 6.2 | 13 | 5.2 | 28 | `radialpayment` | result |
| 7.1 | 15 | 6.1 | 32 | `conditional-law` | result |
| 8.1 | 16 | 6.2 | 33 | `contact-square` | result |
| 8.2 | 17 | 6.3 | 34 | `rem-theta` | result |
| 8.3 | 17 | 6.4 | 35 | `conditional8` | result |
| 9.1 | 18 | 6.5 | 35 | `spectral-reserve` | result |
| 9.2 | 18 | 6.6 | 36 | `rem-radii` | result |
| 9.3 | 19 | 6.7 | 36 | `final-certificate` | result |

| preprint section | thesis location | label |
|---|---|---|
| §2 (p.3) | 3 (p.13) | `sec-variational` |
| §3 (p.6) | 4 (p.18) | `sec-finite` |
| §4 (p.8) | 4.4 (p.21) | `sec-separation` |
| §5 (p.10) | 5.1 (p.25) | `sec-moments` |
| §6 (p.12) | 5.2 (p.27) | `sec-quartic` |
| §7 (p.14) | 6.1 (p.30) | `sec-directional` |
| §8 (p.16) | 6.2 (p.33) | `sec-reserve` |
| §9 (p.18) | 6.3 (p.35) | `sec-conclusion` |
| §10 (p.20) | 8 (p.45) | `sec-further` |
| §A (p.20) | A (p.49) | `sec-appendix` |

### Only in the thesis (new labels)

| thesis number | page | label |
|---|---|---|
| 2.3 | 11 | `bridge` |
| 2.1 | 10 | `blowup-spec` |
| 6.8 | 38 | `weighted-form` |
| 6.10 | 38 | `no-gain` |
| 6.11 | 39 | `padding-note` |
| 7.1 | 41 | `equal-blowup` |
| 1.4 | 2 | `gamma-def` |
| 2.9 | 11 | `dl-bound` |
| 2.6 | 9 | `wg-spectrum` |
| 2.7 | 9 | `cut-eigen` |
| 2.8 | 10 | `quotient-unequal` |
| 2.4 | 9 | `sigmak-def` |
| 2.5 | 9 | `M-def` |
| 2.12 | 12 | `g34` |
| 1.1 | 2 | `hong-bound` |
| 1.3 | 2 | `mohar-bound` |
| 1.2 | 2 | `nik-conj` |
| 1.5 | 2 | `stw-individual` |
| 1.6 | 2 | `stw-hist` |
| 2.10 | 11 | `gammar-def` |
| 2.11 | 11 | `bridge-bound` |
| 2.3 | 8 | `kyfan` |
| 2.1 | 8 | `kn-lambda3` |
| 2.2 | 8 | `k1c5-lambda3` |
| 6.9 | 38 | `weighted-scope` |
| 2.4 | 11 | `bridge-scope` |
| 2.2 | 10 | `unequal-blowup` |

### Missing in thesis aux: none
