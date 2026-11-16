# Cross-references from pivoting-conditioning

- pivoting-conditioning/lec.qmd:24 -> lu (About: solving Ax = b by elimination, assumed known) [ok, linked]
- pivoting-conditioning/lec.qmd:39 -> error-conditioning (notes: same two-source structure, unstable algorithm vs ill-conditioned problem) [ok, linked]
- pivoting-conditioning/lec.qmd:44 -> lu (notes: LU, multipliers, forward/back substitution assumed) [ok, linked]
- pivoting-conditioning/lec.qmd:51 -> error-conditioning (slide "Where we left off": 2x2 two-dye teaser, A = [1 0.95; 0.95 1], b = (1.95, 1.95) -> (1.96, 1.94), x = (1,1) -> (1.2, 0.8), kappa = 39) [ok: same A, b, perturbed b, x and kappa = 39 as error-conditioning slide "Teaser: the condition number of a matrix" and exercise B9]
- pivoting-conditioning/lec.qmd:102 -> least-squares-qr (notes: spectral cytometers with more detectors than dyes = least squares) [ok, linked]
- pivoting-conditioning/lec.qmd:123 -> lu (slide "Recap: elimination as LU": cost 2/3 n^3 factorisation, 2n^2 per right-hand side) [ok, linked]
- pivoting-conditioning/lec.qmd:126 -> lu (notes: one-slide recap only) [ok, linked]
- pivoting-conditioning/lec.qmd:391 -> svd (notes, slide "Matrix norms": singular values explained there) [ok, linked]
- pivoting-conditioning/lec.qmd:425 -> error-conditioning (slide "The condition number of a matrix": kappa = |x f'/f| of a scalar problem) [fixed: "Part 1" -> link to error-conditioning]
- pivoting-conditioning/lec.qmd:504 -> error-conditioning (notes, slide "Rule of thumb": rounding vs measurement error asymmetry, lose log10 kappa digits) [ok, linked]
- pivoting-conditioning/lec.qmd:599 -> least-squares-qr (slide "What to do": more detectors than dyes) [ok, linked]
- pivoting-conditioning/lec.qmd:601 -> svd (slide "What to do": regularisation) [ok, linked]
- pivoting-conditioning/lec.qmd:609 -> least-squares-qr, svd (notes: numerical remedies in Part 4) [ok, linked]
- pivoting-conditioning/lec.qmd:627 -> error-conditioning (notes, checklist: kappa_num, perturb inputs and re-solve) [ok, linked]
- pivoting-conditioning/lec.qmd:649 -> structured-sparse (Summary: next unit) [ok, linked]
- pivoting-conditioning/lec.qmd:657 -> least-squares-qr (notes: normal equations square kappa) [ok, linked]
- pivoting-conditioning/lec.qmd:658 -> identifiability (notes: conditioning returns for fitted parameters) [ok, linked]
- pivoting-conditioning/exercise.qmd:33 -> lu (Background fact 1: elimination, multipliers) [ok, linked]
- pivoting-conditioning/exercise.qmd:984 -> structured-sparse (Where this leads: Cholesky, banded, sparse) [ok, linked]
- pivoting-conditioning/exercise.qmd:985 -> least-squares-qr (Where this leads: unmixing with more detectors than dyes; normal equations square kappa) [fixed: link; added that least-squares-qr teaches the method on calibration curves, not spectra]
- pivoting-conditioning/exercise.qmd:986 -> svd (Where this leads: singular values, regularisation) [ok, linked]
- pivoting-conditioning/exercise.qmd:987 -> newton-systems (Where this leads: each Newton step solves a linear system with the Jacobian) [ok, linked]
- pivoting-conditioning/exercise.qmd:987 -> identifiability (Where this leads: nearly singular matrix = parameters the data cannot pin down) [ok, linked]
- pivoting-conditioning/quiz.qmd:130 (Q10) -> least-squares-qr (more detectors than dyes, "Part 4") [ok, linked]
- pivoting-conditioning/quiz.qmd:130 (Q10) -> svd (regularisation) [ok, linked]

## Running examples defined here

- **Two-dye unmixing matrix** (continued from error-conditioning): $A = [1\ c;\ c\ 1]$,
  rows = detectors at the two emission peaks, columns = dyes, overlap $c$;
  $\kappa_2 = (1+c)/(1-c)$. Standard values: $c = 0.95$ ($\kappa = 39$, peaks about
  10 nm apart for Gaussian spectra of width 30 nm) and $c = 0.5$ ($\kappa = 3$, peaks
  35 nm apart). $\mathbf b = (1.95, 1.95)$, $\mathbf x = (1, 1)$; perturbed
  $\mathbf b = (1.96, 1.94)$ gives $\mathbf x = (1.2, 0.8)$.
- **Six-dye Gaussian panel**: `panel(n, spacing, width = 30)`, peaks at
  `500 + spacing * (0:(n-1))` nm, detectors at the peaks,
  $A_{ij} = \exp(-(p_i - p_j)^2/(2\cdot 30^2))$; true abundances
  $\mathbf x = (5, 1, 3, 2, 4, 1)$. $\kappa$ = 1.64 (60 nm), 28.8 (30 nm), 705 (20 nm),
  $6.4\times10^5$ (10 nm), $1.14\times10^8$ (6 nm), $7.2\times10^8$ (5 nm).
  Natural to reuse in least-squares-qr (more detectors than dyes: put detectors on a
  finer wavelength grid) and svd (regularised unmixing).
- **Noise model for unmixed cells**: 1000 cells with 1 unit of each dye, independent
  Gaussian detector noise sd 0.01 (`set.seed(2)` in R); unmixed sd 0.015 (c = 0.5)
  and 0.142 (c = 0.95).

## Pass 2 (other fixes)

- [fixed: exercise C2, C5, C6 R/Python tabs get a pointer to rel() (C1) and panel() (C5)]
