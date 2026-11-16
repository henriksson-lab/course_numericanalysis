# Cross-references from `fourier-fft` (unit 17)

- fourier-fft/lec.qmd:24 (About this lecture) -> polynomial-interpolation (interpolation with polynomials; previous unit of Part 6) [fixed: linked]
- fourier-fft/lec.qmd:24 (About this lecture) -> splines (interpolation with splines; previous unit) [fixed: linked]
- fourier-fft/lec.qmd:33 (notes, About this lecture) -> sister course quantbio-3MB00A `fourier` (sum of waves, change of basis, masking the spectrum; external, not a unit here) [ok: external, left plain]
- fourier-fft/lec.qmd:52 (Recap: interpolation in a basis) -> polynomial-interpolation (interpolant as a combination of basis functions, Runge phenomenon on equispaced points) [fixed: linked; Runge on equispaced points and basis-function view confirmed in polynomial-interpolation]
- fourier-fft/lec.qmd:150-157 (Trigonometric interpolation in code, bullet + notes) -> least-squares-qr (least-squares fit; normal equations become diagonal for orthogonal columns) [fixed: notes now name Least squares: normal equations vs QR; exercise C3 solution links it]
- fourier-fft/lec.qmd:162 (When the data are not periodic, Check) -> splines (use splines for non-periodic data) [fixed: slide Check now links Splines]
- fourier-fft/lec.qmd:165 (notes, When the data are not periodic) -> convergence-cost (shared logistic OD growth curve K = 1.2, N0 = 0.05, r = 0.8 per hour) [ok: K = 1.2, N0 = 0.05, r = 0.8/h matches convergence-cost and make-fourier-fft-figs.R; notes now name the unit]
- fourier-fft/lec.qmd:247 (notes, O(n^2) against O(n log n)) -> convergence-cost (reading cost slopes off a log-log plot) [ok: log-log cost slopes are in convergence-cost]
- fourier-fft/lec.qmd:390 (notes, Filtering fails at the edges) -> splines (smoothing splines as an alternative smoother) [ok: splines has smoothing splines; notes say 'smoothing splines (Splines)']
- fourier-fft/lec.qmd:422 (Image blur is a convolution) -> svd (regularisation; small singular values amplify noise) [fixed: linked; svd covers Tikhonov/ridge and small singular values]
- fourier-fft/lec.qmd:426 (notes, Image blur) -> svd (Fourier coefficients of the PSF = singular values of the blur matrix) [ok: notes; svd unit]
- fourier-fft/lec.qmd:452 (Summary, Next) -> finite-differences (start of Part 7; differentiation amplifies noise) [fixed: linked; finite-differences is unit 18, first of Part 7, and shows noise amplification on plate-reader data]
- fourier-fft/exercise.qmd:C3(b) and Solution — C3 -> least-squares-qr (cosinor least squares = projection, orthogonal columns) [fixed: linked in C3 solution]
- fourier-fft/exercise.qmd:Solution — C6 ("same bias–noise trade-off as for smoothing splines") -> splines [fixed: linked; splines has 'The bias–noise trade-off' slide]
- fourier-fft/exercise.qmd:1256 (Where this leads) -> finite-differences (differentiation multiplies X_k by 2 pi i k / T) [fixed: linked; note: finite-differences has no Fourier view, but the bullet only claims the shared noise problem and smoothing first, which it does show]
- fourier-fft/exercise.qmd:1257 (Where this leads) -> quadrature (trapezoid rule very accurate for periodic functions over a period) [fixed: linked; note: quadrature shows trapezoid O(h^2) but not the periodic case — reworded as 'the trapezoid rule of that unit is, for periodic functions, far more accurate than its usual O(h^2)']
- fourier-fft/exercise.qmd:1258 (Where this leads) -> diffusion-explicit (Fourier modes decay at own rates; stability set by the Nyquist mode) [fixed: linked; diffusion-explicit has the Fourier sine-series solution, modes decaying at own rates, and the zigzag instability]
- fourier-fft/exercise.qmd:1259 (Where this leads) -> svd (deconvolution as ill-conditioned problem) [fixed: linked]
- fourier-fft/quiz.qmd:Q12 (Answer) -> svd (small singular values, ridge = Wiener filter) [fixed: linked in Q12 answer]

## Running examples defined here

- **Circadian clock-controlled gene** (lec, exercise C3, figures). 24 samples, every 2 h over 48 h (t = 0, 2, ..., 46),
  log2 expression, generated once with set.seed(1) as 8 + 1.5 cos(2 pi (t - 6)/24) + N(0, 0.3^2), rounded to 2 decimals:
  x = (7.81, 8.81, 9.05, 9.98, 9.40, 8.50, 8.15, 7.47, 6.87, 6.41, 7.15, 7.37,
       7.81, 8.09, 9.64, 9.49, 9.29, 9.03, 8.25, 7.43, 6.98, 6.73, 6.72, 6.65).
  Key numbers: mean 8.045, 24-h amplitude 1.4645, peak 6.45 h; cosinor coefficients (8.0450, -0.1731, 1.4542).
  Could be reused in `nonlinear-least-squares` (fit period as a parameter) or `identifiability`.
- **Ultradian alias**: cos(2 pi t / 4.8) sampled every 4 h is identical to a 24-h cosine.
- **PER2::LUC-like bioluminescence trace**: t = 0, 0.5, ..., 119.5 h (n = 240),
  truth = 20 + 60 exp(-t/40) + 30 exp(-t/80) cos(2 pi (t - 4)/24), noise sd 6 (R set.seed(3)).
  Could be reused in `finite-differences` (derivative of a noisy trace) or `splines`.
