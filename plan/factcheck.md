# Facts to verify against sources (teacher-only)

Claims that workers wrote from memory or without reading the source. Needs web access, so it is left for when the course lead is present (web lookups may need approval). Check each, then fix or soften the text in the named unit.

- [ ] pivoting-conditioning: Cytek "complexity index" = condition number of the reference-spectra matrix — stated from memory, verify or soften.
- [ ] svd: tool names EPIC, quanTIseq, MuSiC (slide) and CIBERSORT (notes) as constrained-LS deconvolution — from memory, verify. Only 4/27 slides screenshotted. scipy 1.15.2 nnls bug noted in C6 (may go stale). Install: nnls (R).
- [ ] finite-differences: optim default BFGS gradient = central diff with ndeps=1e-3 — verify. Install: numDeriv.
- [ ] splines: tradeSeq / Monocle use spline trends — from memory, verify.
- [ ] monte-carlo: from memory — RANDU history, R 3.6 sample() change, ChaCha12 as StdRng, Bellman 'curse of dimensionality', bedtools shuffle / GAT / regioneR. Slides 13–18 figures not seen in browser.
- [ ] runge-kutta: Python C2/C3/C5/C6 tabs import C1/C4 functions from local files (not self-contained) — make self-contained in audit. Repressilator params α=216, α0=0.216, β=5, n=2 unverified vs Elowitz & Leibler. ode_solvers Dopri5 dense-output end-point bug noted.
- [ ] fourier-fft: from memory — 4.8-h ultradian example, p53 ~5.5 h, Hes1 ~2 h periods, 'sample every 2 h' guideline. R fft prime-length slowness claim measured here.
- [ ] constrained-lp: from memory — iML1515 2712 reactions / 1877 metabolites, Recon3D ~13 000 reactions. Install: lpSolve, numDeriv.
- [ ] quadrature: from memory — bioequivalence 80–125 %, regulators accept sampled AUC, WinNonlin/PKNCA rules. Romberg recap assumes Richardson (4D(h/2)−D(h))/3 form — check vs finite-differences in pass 2.
- [ ] gradient-methods: slides not screenshotted; densest 'What optim, minimize and argmin do', 'When Newton fails'. From memory: optim ndeps=1e-3 central; scipy forward diff step √ε·max(1,|x|); BFGS gtol 1e-5; McKinnon 1998. rust/examples/gradient_methods_shared/defs.rs included via include! (new pattern).
- [ ] identifiability: pk(t, ka, ke, V, D=100) vs NLS unit's bateman() naming — reconcile in pass 2. From memory: Raue 2009, Gutenkunst 2007, DAISY, StructuralIdentifiability.jl. Slides not screenshotted.
- [ ] steady-state-2d: [B1] 'mechanics not nutrients' claim from contentplan, not read — verify. Grid conventions: x_i = i h, µm/s/µM, ghost points, k = i + (j−1)nx — reconcile with diffusion-implicit in pass 2.
- [ ] diffusion-implicit: R/Python solution tabs rely on Background-defined helpers (thomas, be_step, cn_step, exact) — decide self-contained convention in audit. h = 10 W/(m²K) from memory. Read [M2], [H1] full text (agent fetched before no-network rule); [M1], [H2] headline only. Embryo 1D model is our own.
- [ ] fem-outlook: notes assume content of steady-state-2d, reaction-diffusion, iterative-solvers, advection-diffusion lecture (T-sensor 100 µm, 1000 µm/s) — verify in pass 2. Hydrogel numbers own round values.
- [ ] advection-diffusion: cell-centred grid across channel (vs nodes elsewhere?) — reconcile with steady-state-2d. Device numbers own round values; Kamholz & Yager 2001 cited cautiously. Some code blocks 13–14 lines.
- [ ] reaction-diffusion: time in minutes vs seconds in diffusion-explicit — reconcile. 2D uses periodic FFT diffusion, not Kronecker Laplacian. From memory: Bicoid λ≈100 µm, D FRAP 0.3 / FCS 7 µm²/s; Fgf8 D≈50, λ≈200 µm; Gregor 2007, Abu-Arish 2010, Berezhkovskii 2010 (notes only). Slides not screenshotted.

## Part 10 literature

- [ ] Check every number attributed to [G1]–[M4] in units 31–35 against the papers listed in contentplan.md (only [M2] and [H1] were read in full).
- [ ] sde: from memory — water viscosity 1.0 mPa s at 20 °C and k_BT = 4.05e-21 J (Stokes–Einstein bead D = 0.43 µm²/s); optical-trap stiffness ~0.01 pN/nm giving theta ≈ 1/ms; OU trait evolution (Hansen 1997; Butler & King 2004); Lewontin & Cohen 1969 (random-environment growth); packages sde, Sim.DiffProc (R), sdeint (Python). Slides with figures not screenshotted (headless Chrome did not load images; r-stretch should fit them).
