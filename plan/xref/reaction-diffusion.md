# Cross-references from reaction-diffusion

- reaction-diffusion/lec.qmd:33 (notes, About) -> structured-sparse, stiffness, eigen-algorithms, phase-plane, fourier-fft (Parts 3, 8, 5, 6 by part number only) [fixed: linked units; added sentence: earlier Part 10 units used seconds, 100 µm²/min = 1.7 µm²/s]
- reaction-diffusion/lec.qmd:51 -> structured-sparse (morphogen in a row of cells, thomas(), decay length 10 cells) [ok: D = 1/min, k = 0.01/min, source cell 1, 10 cells verified; linked]
- reaction-diffusion/lec.qmd:58 -> diffusion-explicit (grid x_i = i dx, (1,-2,1) stencil, ghost points, FTCS limit r <= 1/2) [ok; linked]
- reaction-diffusion/lec.qmd:58 -> diffusion-implicit (backward Euler, one Thomas solve per step) [ok; linked]
- reaction-diffusion/lec.qmd:123 (notes) -> diffusion-implicit (same lo/di/up layout; 1 + dt * diagonals) [fixed: diffusion-implicit does not use lo/di/up names; reworded to "same thomas(a, b, c, d) call"]
- reaction-diffusion/lec.qmd:123 (notes) -> structured-sparse (thomas unchanged; diagonal dominance, no pivoting) [ok; linked]
- reaction-diffusion/lec.qmd:134 (notes) -> structured-sparse (residual / mass-balance checks) [fixed: linked; named the two checks]
- reaction-diffusion/lec.qmd:158 -> stiffness (stiff diffusion, treat implicitly) [ok; linked]
- reaction-diffusion/lec.qmd:163 (notes) -> diffusion-implicit (backward Euler with k on the diagonal); newton-systems (Newton per step for fully implicit nonlinear reaction) [ok; linked]
- reaction-diffusion/lec.qmd:188 (notes) -> diffusion-implicit (same Thomas call as backward Euler) [ok; linked]
- reaction-diffusion/lec.qmd:223 (notes) -> advection-diffusion (splitting also used for advection) [ok; linked]
- reaction-diffusion/lec.qmd:259 (notes) -> diffusion-implicit (temperature step inside the embryo, 20/27 C streams [M1, M2]) [fixed: diffusion-implicit finds a ~130 µm wide step, not "sharp"; reworded]
- reaction-diffusion/lec.qmd:269 (notes) -> least-squares-qr (Part 4: linear least squares), nonlinear-least-squares (Part 9: fit background as parameter) [ok; linked]
- reaction-diffusion/lec.qmd:279 -> phase-plane (Jacobian, trace/determinant stability) [ok; linked]
- reaction-diffusion/lec.qmd:306 (notes) -> diffusion-explicit (Fourier modes of the diffusion equation); phase-plane (2x2 eigenvalue test) [ok; linked]
- reaction-diffusion/lec.qmd:323 (notes) -> phase-plane (trace-determinant classification) [ok; linked]
- reaction-diffusion/lec.qmd:345 (notes) -> steady-state-2d (five-point Laplacian; sparse direct solve with walls) [fixed: added one sentence why periodic FFT instead of the sparse Kronecker Laplacian; linked]
- reaction-diffusion/lec.qmd:345 (notes) -> fourier-fft (FFT diagonalises the periodic Laplacian; plain "fft") [fixed: added one sentence why periodic FFT instead of the sparse Kronecker Laplacian; linked]
- reaction-diffusion/lec.qmd:379 -> euler (Euler stability |1 + h lambda| <= 1) [ok; linked]
- reaction-diffusion/lec.qmd:401 (notes) -> diffusion-explicit (checkerboard = unstable explicit diffusion step) [ok; linked]
- reaction-diffusion/lec.qmd:"What you can now do" slide (l. 418-424) -> Parts 1-10 by part number (course wrap-up) [fixed: rewritten to match sidebar parts (SVD in Part 4, Monte Carlo, stochastic simulation, constraints, 2D steady states, flow added)]
- reaction-diffusion/lec.qmd:428 (notes) -> fem-outlook, iterative-solvers (Part 11 further reading) [ok; linked]
- reaction-diffusion/exercise.qmd:37, 44, 84 -> structured-sparse (thomas()) [ok: thomas defined on page; linked]
- reaction-diffusion/exercise.qmd:229 (B5(b)) -> phase-plane (type of fixed point) [ok; linked]
- reaction-diffusion/exercise.qmd:1066 (C5 solution) -> nonlinear-least-squares (fit c0 exp(-x/lambda) + b) [ok; linked]
- reaction-diffusion/exercise.qmd:Where this leads -> fem-outlook, iterative-solvers, nonlinear-least-squares, identifiability, diffusion-implicit [fixed: linked]

## Grid conventions used here (for reconciliation with units 31-34)

- Nodes x_i = i dx, i = 0..N (as in diffusion-explicit and the diffusion-implicit figure script); ghost points for no-flux and flux boundaries: c_{-1} = c_1 + 2 dx j / D (influx j at x = 0), c_{N+1} = c_{N-1} (wall).
- The operator -D d2/dx2 is stored as three diagonals lo, di, up (positive operator K_D): di = 2r, off-diagonals -r, and -2r in row 0 (above) and row N (below), r = D/dx^2; influx right-hand side 2 j / dx in row 0. Thomas convention as in structured-sparse (a[1], c[n] unused).
- Units: µm and **minutes** (to match the structured-sparse row of cells), not seconds as in diffusion-explicit.
- steady-state-2d and advection-diffusion did not exist when this unit was written; the 2D part here uses a periodic grid and FFT, not the Kronecker-product sparse Laplacian.

## Running examples defined here

- **Bicoid-like SDD gradient** (continuum version of the structured-sparse row of cells, cell = 10 µm): D = 100 µm^2/min, k = 0.01 /min, L = 500 µm, influx j = D/lambda = 1 at x = 0, wall at x = L; lambda = 100 µm. Exact steady state c = (j lambda / D) cosh((L-x)/lambda)/sinh(L/lambda), c(0) = 1.0001. Grid N = 50 (dx = 10 µm) max error 0.001247396; accumulation times 50, 100, 150 min at x = 0, 100, 200 µm.
- **Temperature step (assumed)**: k = 0.01 cold / 0.02 warm (k doubles at 27 C), D and j unchanged; threshold c = 0.1 at 230.7 (cold), 138.3 (warm), 220.2 (posterior warm), 138.9 µm (anterior warm).
- **Schnakenberg** a = 0.1, b = 0.9, d = 20 (D_u = 1): uniform state (1, 0.9), J = [0.8 1; -1.8 -1]; unstable band q in [0.272, 0.822], fastest q = 0.487 (wavelength 12.9), d_c = 8.5676. 2D: periodic 80 x 80, dx = 1, IMEX dt = 0.1 to t = 200 -> 47 spots. [ok; linked]
- **lcg(n, seed)**: Park-Miller generator s <- (16807 s) mod (2^31 - 1), returns s / (2^31 - 1); identical numbers in R, Python, Rust.
