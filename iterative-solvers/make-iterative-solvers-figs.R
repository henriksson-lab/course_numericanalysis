# Figures for the iterative-solvers/ unit ("Iterative linear solvers in depth").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd iterative-solvers && Rscript make-iterative-solvers-figs.R
#
# Base R graphics + Matrix only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Takes about five minutes (the 3D Cholesky factorisations
# and the slow Jacobi / Gauss-Seidel runs dominate). Prints the numbers quoted
# on the slides.
#
# Running example: oxygen in a cell-laden gel in a culture chamber, square
# cross-section L = 1 mm, oxygen supplied through all four walls (c = 1, scaled
# by the ambient concentration), uniform consumption. Unknown: the oxygen
# deficit u = 1 - c, which solves -Laplace(u) = Phi with u = 0 on the walls,
# Phi = q L^2 / (D c0) = 12. On an m x m interior grid, h = 1/(m+1), multiplied
# through by h^2: A u = Phi h^2 * 1, A = 2D five-point Laplacian (4 on the
# diagonal, -1 for the four neighbours). Tolerance: relative residual 1e-6.
#
# Makes:
#   fig01-chamber.svg     the oxygen field and the hypoxic core (m = 127)
#   fig02-fillin.svg      non-zeros of the Cholesky factor vs unknowns, 2D and 3D
#   fig03-smooth.svg      Gauss-Seidel error after 0, 10, 200 sweeps: it becomes smooth
#   fig04-residuals.svg   residual vs iteration, m = 63: Jacobi, GS, CG, IC-CG
#   fig05-counts.svg      iterations vs grid size, log-log, with multigrid
#   fig06-gmres.svg       perfused chamber (advection-diffusion): field, and
#                         residuals of CG, GMRES and GMRES(20)
#   fig07-multigrid.svg   multigrid in one picture: smooth on the fine grid,
#                         correct on coarse grids (V-cycle)

library(Matrix)
oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999", black = "#000000")
par_big <- function(...) par(cex.lab = 1.5, cex.axis = 1.35, cex.main = 1.5, ...)
seq_pal <- colorRampPalette(c("#2D004B", "#542788", "#8073AC", "#B2ABD2",
                              "#FDB863", "#E08214", "#FFF5E0"))(100)  # dark = low O2

## ---- matrices -----------------------------------------------------------------

lap1 <- function(m) bandSparse(m, k = -1:1,
  diagonals = list(rep(-1, m - 1), rep(2, m), rep(-1, m - 1)))
lap2 <- function(m) { I <- Diagonal(m); T <- lap1(m); kronecker(I, T) + kronecker(T, I) }
lap3 <- function(m) { I <- Diagonal(m); T <- lap1(m)
  kronecker(kronecker(I, I), T) + kronecker(kronecker(I, T), I) + kronecker(kronecker(T, I), I) }
# perfused chamber: flow in +x with Peclet number Pe, upwind differences
adv2 <- function(m, Pe) { h <- 1 / (m + 1)
  Tx <- lap1(m) + Pe * h * bandSparse(m, k = c(0, -1), diagonals = list(rep(1, m), rep(-1, m - 1)))
  kronecker(Diagonal(m), Tx) + kronecker(lap1(m), Diagonal(m)) }
# variable diffusion coefficient Dfun(x, y), evaluated on the cell faces
lapD <- function(m, Dfun) {
  h <- 1 / (m + 1); I <- Diagonal(m)
  G1 <- bandSparse(m + 1, m, k = c(0, -1), diagonals = list(rep(1, m), rep(-1, m)))
  Gx <- kronecker(I, G1); Gy <- kronecker(G1, I)
  fx <- expand.grid(x = (1:(m + 1) - 0.5) * h, y = (1:m) * h)        # x-faces
  fy <- expand.grid(x = (1:m) * h, y = (1:(m + 1) - 0.5) * h)        # y-faces
  t(Gx) %*% Diagonal(x = Dfun(fx$x, fx$y)) %*% Gx +
    t(Gy) %*% Diagonal(x = Dfun(fy$x, fy$y)) %*% Gy
}
aggregate_D <- function(x, y) ifelse((x - 0.5)^2 + (y - 0.5)^2 < 0.3^2, 0.1, 1)

## ---- solvers (same code as on the slides / in the exercise) --------------------

relres <- function(A, b, x) sqrt(sum((b - A %*% x)^2)) / sqrt(sum(b^2))
jacobi <- function(A, b, tol = 1e-6, maxit = 1e6, trace = FALSE) {
  x <- numeric(length(b)); d <- diag(A); nb <- sqrt(sum(b^2)); res <- numeric(0)
  for (it in 0:maxit) {
    r <- b - (A %*% x)[, 1]; rel <- sqrt(sum(r^2)) / nb
    if (trace) res <- c(res, rel)
    if (rel < tol) break
    x <- x + r / d
  }
  list(x = x, it = it, res = res)
}
gauss_seidel <- function(A, b, tol = 1e-6, maxit = 1e6, trace = FALSE) {
  x <- numeric(length(b)); M <- tril(A); nb <- sqrt(sum(b^2)); res <- numeric(0)
  for (it in 0:maxit) {
    r <- b - (A %*% x)[, 1]; rel <- sqrt(sum(r^2)) / nb
    if (trace) res <- c(res, rel)
    if (rel < tol) break
    x <- x + as.vector(solve(M, r))
  }
  list(x = x, it = it, res = res)
}
pcg <- function(A, b, Minv = function(r) r, tol = 1e-6, maxit = 1e5) {
  x <- numeric(length(b)); r <- b; z <- Minv(r); p <- z; rz <- sum(r * z)
  nb <- sqrt(sum(b^2)); res <- 1
  for (it in 1:maxit) {
    Ap <- (A %*% p)[, 1]; a <- rz / sum(p * Ap)
    x <- x + a * p; r <- r - a * Ap
    res <- c(res, sqrt(sum(r^2)) / nb)
    if (res[it + 1] < tol) break
    z <- Minv(r); rz_new <- sum(r * z); p <- z + (rz_new / rz) * p; rz <- rz_new
  }
  list(x = x, it = it, res = res)
}
# CG that also records the TRUE residual (for the nonsymmetric failure)
cg_true <- function(A, b, maxit = 300) {
  x <- numeric(length(b)); r <- b; p <- r; rr <- sum(r^2); res <- 1
  for (it in 1:maxit) {
    Ap <- (A %*% p)[, 1]; a <- rr / sum(p * Ap); x <- x + a * p; r <- r - a * Ap
    res <- c(res, relres(A, b, x)); rn <- sum(r^2); p <- r + (rn / rr) * p; rr <- rn
  }
  res
}
# IC(0) for a five-point matrix in natural ordering: L has the pattern of tril(A)
ic0 <- function(A, m) {
  n <- nrow(A); d <- numeric(n); lw <- numeric(n); ls <- numeric(n)
  dA <- diag(A); w <- c(0, diag(A[-1, -n])); s <- c(rep(0, m), diag(A[-(1:m), -((n - m + 1):n)]))
  for (k in 1:n) {
    if (k > 1) lw[k] <- w[k] / d[k - 1]
    if (k > m) ls[k] <- s[k] / d[k - m]
    d[k] <- sqrt(dA[k] - lw[k]^2 - ls[k]^2)
  }
  bandSparse(n, k = c(0, -1, -m), diagonals = list(d, lw[-1], ls[-(1:m)]))
}
ic_prec <- function(A, m) { L <- ic0(A, m); Lt <- t(L)
  function(r) as.vector(solve(Lt, solve(L, r))) }
# full GMRES (no restart), residual history from the small least-squares problem
gmres <- function(A, b, k = 200, tol = 1e-6, x0 = numeric(length(b))) {
  r0 <- b - (A %*% x0)[, 1]; n <- length(b)
  Q <- matrix(0, n, k + 1); H <- matrix(0, k + 1, k)
  beta <- sqrt(sum(r0^2)); Q[, 1] <- r0 / beta; res <- numeric(0)
  for (j in 1:k) {
    w <- (A %*% Q[, j])[, 1]
    for (i in 1:j) { H[i, j] <- sum(Q[, i] * w); w <- w - H[i, j] * Q[, i] }
    H[j + 1, j] <- sqrt(sum(w^2)); Q[, j + 1] <- w / H[j + 1, j]
    e1 <- c(beta, rep(0, j)); Hj <- H[1:(j + 1), 1:j, drop = FALSE]
    y <- qr.solve(Hj, e1)
    res <- c(res, sqrt(sum((e1 - Hj %*% y)^2)) / sqrt(sum(b^2)))
    if (tail(res, 1) < tol) break
  }
  list(x = x0 + as.vector(Q[, 1:j] %*% y), it = j, res = res)
}
gmres_restart <- function(A, b, k = 20, tol = 1e-6, maxcycle = 100) {
  x <- numeric(length(b)); res <- 1
  for (cyc in 1:maxcycle) {
    g <- gmres(A, b, k = k, tol = tol, x0 = x); x <- g$x
    res <- c(res, g$res)                  # g$res is relative to ||b||
    if (tail(res, 1) < tol) break
  }
  list(x = x, res = res)
}
# geometric multigrid V(2,2) with Gauss-Seidel smoothing, Galerkin coarse grids
prol1 <- function(mc) {
  mf <- 2 * mc + 1
  sparseMatrix(i = c(2 * (1:mc) - 1, 2 * (1:mc), 2 * (1:mc) + 1),
               j = rep(1:mc, 3), x = rep(c(0.5, 1, 0.5), each = mc), dims = c(mf, mc))
}
mg_setup <- function(A, m) {
  lev <- list()
  while (m > 3) {
    mc <- (m - 1) / 2; P <- kronecker(prol1(mc), prol1(mc))
    lev[[length(lev) + 1]] <- list(A = A, P = P, M = tril(A))
    A <- t(P) %*% A %*% P / 4; m <- mc
  }
  lev[[length(lev) + 1]] <- list(A = A)
  lev
}
vcycle <- function(lev, l, b, x) {
  A <- lev[[l]]$A
  if (l == length(lev)) return(as.vector(solve(A, b)))
  M <- lev[[l]]$M; P <- lev[[l]]$P
  for (s in 1:2) x <- x + as.vector(solve(M, b - (A %*% x)[, 1]))
  r <- b - (A %*% x)[, 1]
  x <- x + as.vector(P %*% vcycle(lev, l + 1, as.vector(t(P) %*% r) / 4, numeric(ncol(P))))
  for (s in 1:2) x <- x + as.vector(solve(M, b - (A %*% x)[, 1]))
  x
}
multigrid <- function(A, b, m, tol = 1e-6, maxit = 100) {
  lev <- mg_setup(A, m); x <- numeric(length(b))
  for (it in 1:maxit) { x <- vcycle(lev, 1, b, x); if (relres(A, b, x) < tol) break }
  it
}

Phi <- 12
rhs <- function(m) rep(Phi / (m + 1)^2, m^2)

## ---- fig01: the chamber -----------------------------------------------------------

m <- 127; A <- lap2(m); u <- solve(A, rhs(m)); cfield <- 1 - matrix(as.vector(u), m, m)
g <- (1:m) / (m + 1)
cat("centre O2 (m = 127):", cfield[64, 64], " fraction below 0.25:", mean(cfield < 0.25), "\n")
svg("fig01-chamber.svg", width = 11, height = 5.6)
layout(matrix(1:2, 1), widths = c(1, 1.15))
par_big(mar = c(4.5, 5, 2.5, 1))
plot(NA, xlim = c(0, 1), ylim = c(0, 1), asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = "Cell-laden gel, 1 mm x 1 mm")
rect(0, 0, 1, 1, col = "#F3EEE2", border = oi["blue"], lwd = 8)
set.seed(3); points(runif(250), runif(250), pch = 19, col = oi["orange"], cex = 0.8)
text(0.5, 1.07, expression(O[2] ~ "in through the walls:" ~ c == 1), cex = 1.3, xpd = NA)
text(0.5, 0.5, "every cell\nconsumes O2", cex = 1.4, font = 2)
par_big(mar = c(4.5, 5, 2.5, 6))
image(g, g, cfield, col = seq_pal, zlim = c(0, 1), asp = 1, xlab = "x (mm)", ylab = "y (mm)",
      main = "Oxygen c = 1 - u  (16 129 unknowns)")
contour(g, g, cfield, levels = 0.25, add = TRUE, lwd = 3, col = "white", labcex = 1.2,
        drawlabels = FALSE)
text(0.5, 0.5, "hypoxic\nc < 0.25", col = "white", cex = 1.3, font = 2)
par(new = TRUE, mar = c(4.5, 0, 2.5, 1.5), fig = c(0.93, 0.97, 0.18, 0.86))
image(1, seq(0, 1, length.out = 100), t(seq(0, 1, length.out = 100)), col = seq_pal,
      axes = FALSE, xlab = "", ylab = "")
axis(4, cex.axis = 1.2, las = 1)
dev.off()

## ---- fig02: fill-in in 2D and 3D ---------------------------------------------------

fill <- data.frame()
for (m in c(30, 50, 100, 200, 400, 800)) {
  A <- lap2(m); t <- system.time(L <- Cholesky(A, perm = TRUE, LDL = FALSE))[3]
  fill <- rbind(fill, data.frame(dim = 2, m = m, n = m^2, nnzA = nnzero(A),
                                 nnzL = nnzero(as(L, "Matrix")), time = t))
}
for (m in c(10, 15, 20, 30, 40, 50)) {
  A <- lap3(m); t <- system.time(L <- Cholesky(A, perm = TRUE, LDL = FALSE))[3]
  fill <- rbind(fill, data.frame(dim = 3, m = m, n = m^3, nnzA = nnzero(A),
                                 nnzL = nnzero(as(L, "Matrix")), time = t))
}
print(fill)
svg("fig02-fillin.svg", width = 11, height = 5.4)
par_big(mar = c(5, 5.5, 1, 1))
plot(NA, log = "xy", xlim = c(800, 1e6), ylim = c(4e3, 5e7), xaxt = "n", yaxt = "n",
     xlab = "unknowns n", ylab = "stored non-zeros")
axis(1, at = 10^(3:6), labels = c("1e3", "1e4", "1e5", "1e6"))
axis(2, at = 10^(4:7), labels = c("1e4", "1e5", "1e6", "1e7"), las = 1)
for (d in 2:3) {
  s <- fill[fill$dim == d, ]; col <- if (d == 2) oi["blue"] else oi["red"]
  lines(s$n, s$nnzL, type = "b", pch = 19, lwd = 3, col = col, cex = 1.3)
  lines(s$n, s$nnzA, type = "b", pch = 1, lwd = 2, lty = 2, col = col, cex = 1.3)
}
legend("topleft", bty = "n", cex = 1.3, lwd = c(3, 2, 3, 2), lty = c(1, 2, 1, 2),
       pch = c(19, 1, 19, 1), col = oi[c("red", "red", "blue", "blue")],
       legend = c("3D: Cholesky factor L", "3D: matrix A (7-point)",
                  "2D: Cholesky factor L", "2D: matrix A (5-point)"))
dev.off()

## ---- fig03: the Gauss-Seidel error becomes smooth ----------------------------------
# (plain Jacobi also keeps a checkerboard mode, eigenvalue -cos(pi h); see notes)

m <- 31; A <- lap2(m); b <- rhs(m); xs <- as.vector(solve(A, b)); g <- (1:m) / (m + 1)
set.seed(1); x <- xs + rnorm(m^2) * 0.02; M <- tril(A)
errs <- list()
for (it in 0:200) {
  if (it %in% c(0, 10, 200)) errs[[as.character(it)]] <- matrix(x - xs, m, m)
  x <- x + as.vector(solve(M, b - (A %*% x)[, 1]))
}
svg("fig03-smooth.svg", width = 13, height = 4.6)
par(mfrow = c(1, 3)); par_big(mar = c(2, 2, 3.5, 1))
div <- colorRampPalette(c(oi["blue"], "white", oi["red"]))(101)
for (nm in names(errs)) {
  e <- errs[[nm]]; z <- max(abs(e))
  image(g, g, e, col = div, zlim = c(-z, z), axes = FALSE, asp = 1, xlab = "", ylab = "",
        main = sprintf("after %s sweeps: max %.1e", nm, z))
  box()
}
dev.off()

## ---- fig04: residual vs iteration, m = 64 ------------------------------------------

m <- 63; A <- lap2(m); b <- rhs(m)
rj <- jacobi(A, b, trace = TRUE); rg <- gauss_seidel(A, b, trace = TRUE)
rc <- pcg(A, b); ri <- pcg(A, b, ic_prec(A, m)); rd <- pcg(A, b, function(r) r / diag(A))
cat("m = 63: Jacobi", rj$it, " GS", rg$it, " CG", rc$it, " Jacobi-PCG", rd$it, " IC-CG", ri$it, "\n")
svg("fig04-residuals.svg", width = 11, height = 5.2)
par_big(mar = c(5, 5.5, 1, 1))
plot(NA, log = "xy", xlim = c(1, 2e4), ylim = c(1e-6, 1.5), yaxt = "n",
     xlab = "iteration (log scale)", ylab = "relative residual")
axis(2, at = 10^(-6:0), labels = c("1e-6", "1e-5", "1e-4", "1e-3", "0.01", "0.1", "1"), las = 1)
k <- function(r) seq_along(r)
lines(k(rj$res), rj$res, lwd = 3, col = oi["orange"])
lines(k(rg$res), rg$res, lwd = 3, col = oi["purple"])
lines(k(rc$res), rc$res, lwd = 3, col = oi["blue"])
lines(k(ri$res), ri$res, lwd = 3, col = oi["green"])
legend("bottomleft", bty = "n", cex = 1.3, lwd = 3, col = oi[c("orange", "purple", "blue", "green")],
       legend = c(sprintf("Jacobi (%d)", rj$it), sprintf("Gauss-Seidel (%d)", rg$it),
                  sprintf("CG (%d)", rc$it), sprintf("CG + incomplete Cholesky (%d)", ri$it)))
dev.off()

## ---- fig05: iterations vs grid size --------------------------------------------------

ms <- c(15, 31, 63, 127, 255, 511)
cnt <- data.frame(m = ms, jac = NA, gs = NA, cg = NA, ic = NA, mg = NA)
for (i in seq_along(ms)) {
  m <- ms[i]; A <- lap2(m); b <- rhs(m)
  if (m <= 127) { cnt$jac[i] <- jacobi(A, b)$it; cnt$gs[i] <- gauss_seidel(A, b)$it }
  cnt$cg[i] <- pcg(A, b)$it; cnt$ic[i] <- pcg(A, b, ic_prec(A, m))$it
  cnt$mg[i] <- multigrid(A, b, m)
}
cnt$jac_pred <- ceiling(log(1e-6) / log(cos(pi / (cnt$m + 1))))
print(cnt)
svg("fig05-counts.svg", width = 11, height = 5.4)
par_big(mar = c(5, 5.5, 1, 1))
h <- 1 / (ms + 1)
plot(NA, log = "xy", xlim = c(15, 511), ylim = c(3, 3e5), xaxt = "n", yaxt = "n",
     xlab = "grid points per side m  (unknowns m^2)", ylab = "iterations to 1e-6")
axis(1, at = ms); axis(2, at = 10^(0:5), labels = c("1", "10", "100", "1e3", "1e4", "1e5"), las = 1)
lines(ms, cnt$jac_pred, lty = 3, lwd = 2, col = oi["orange"])
lines(ms, cnt$jac, type = "b", pch = 19, lwd = 3, col = oi["orange"], cex = 1.3)
lines(ms, cnt$gs, type = "b", pch = 15, lwd = 3, col = oi["purple"], cex = 1.3)
lines(ms, cnt$cg, type = "b", pch = 17, lwd = 3, col = oi["blue"], cex = 1.3)
lines(ms, cnt$ic, type = "b", pch = 18, lwd = 3, col = oi["green"], cex = 1.6)
lines(ms, cnt$mg, type = "b", pch = 19, lwd = 3, col = oi["black"], cex = 1.3)
text(255, 6e4, "Jacobi: ~ 1/h^2\n(dotted: from cos(pi h))", col = oi["orange"], cex = 1.2)
text(400, 900, "CG: ~ 1/h", col = oi["blue"], cex = 1.3)
text(400, 7, "multigrid: constant", cex = 1.3)
legend("topleft", bty = "n", cex = 1.2, lwd = 3, pch = c(19, 15, 17, 18, 19),
       col = oi[c("orange", "purple", "blue", "green", "black")],
       legend = c("Jacobi", "Gauss-Seidel", "CG", "CG + incomplete Cholesky", "multigrid V-cycles"))
dev.off()

## preconditioning experiment: uniform and aggregate chambers, m = 64 and 128
for (m in c(63, 127)) {
  for (case in c("uniform", "aggregate")) {
    A <- if (case == "uniform") lap2(m) else lapD(m, aggregate_D); b <- rhs(m)
    cat(case, m, " CG", pcg(A, b)$it, " Jacobi-PCG", pcg(A, b, function(r) r / diag(A))$it,
        " IC-CG", pcg(A, b, ic_prec(A, m))$it, "\n")
  }
}

## ---- fig06: perfused chamber, GMRES ---------------------------------------------------

m <- 63; Pe <- 50; A <- adv2(m, Pe); b <- rhs(m); g <- (1:m) / (m + 1)
cat("symmetric?", isSymmetric(A), "\n")
rcg <- cg_true(A, b, 300); rgm <- gmres(A, b, k = 300); rgr <- gmres_restart(A, b, k = 20)
cat("Pe = 50: CG min res", min(rcg), " last", tail(rcg, 1), "  GMRES its", rgm$it,
    "  GMRES(20) its", length(rgr$res) - 1, "\n")
cat("GMRES on the symmetric chamber (m = 63):", gmres(lap2(m), b, k = 300)$it,
    "  CG:", pcg(lap2(m), b)$it, "\n")
cfield <- 1 - matrix(as.vector(solve(A, b)), m, m)
svg("fig06-gmres.svg", width = 13, height = 5.4)
layout(matrix(1:2, 1), widths = c(1, 1.25))
par_big(mar = c(4.5, 4.5, 3, 1))
image(g, g, cfield, col = seq_pal, zlim = c(0, 1), asp = 1, xlab = "x (mm)", ylab = "y (mm)",
      main = "Perfused chamber, flow ->, Pe = 50")
contour(g, g, cfield, levels = 0.25, add = TRUE, lwd = 3, col = "white", drawlabels = FALSE)
arrows(0.1, 0.5, 0.35, 0.5, lwd = 4, col = "white", length = 0.15)
par_big(mar = c(4.5, 5.5, 3, 1))
plot(NA, xlim = c(0, 300), ylim = c(1e-6, 1e4), log = "y", yaxt = "n",
     xlab = "iteration (matrix-vector products)", ylab = "relative residual", main = "Nonsymmetric A")
axis(2, at = 10^c(-6, -4, -2, 0, 2, 4), labels = c("1e-6", "1e-4", "0.01", "1", "100", "1e4"), las = 1)
lines(seq_along(rcg) - 1, rcg, lwd = 3, col = oi["red"])
lines(seq_along(rgr$res) - 1, rgr$res, lwd = 3, col = oi["sky"])
lines(0:rgm$it, c(1, rgm$res), lwd = 3, col = oi["blue"])
legend("topright", bty = "n", cex = 1.25, lwd = 3, col = oi[c("red", "blue", "sky")],
       legend = c("CG (wrong method)", sprintf("GMRES (%d)", rgm$it),
                  sprintf("GMRES(20), restarted (%d)", length(rgr$res) - 1)))
dev.off()

## ---- fig07: multigrid in one picture -------------------------------------------------

svg("fig07-multigrid.svg", width = 13, height = 5.2)
layout(matrix(1:2, 1), widths = c(1, 1.2))
# left: 1D error before / after smoothing, and on the coarse grid
par_big(mar = c(4.5, 4.5, 3, 1))
m1 <- 63; x1 <- (1:m1) / (m1 + 1); set.seed(2)
e0 <- sin(pi * x1) + 0.5 * rnorm(m1) * sin(pi * x1)
T1 <- as.matrix(lap1(m1)); e <- e0
for (s in 1:5) e <- e - (2 / 3) * (T1 %*% e)[, 1] / 2         # weighted Jacobi on A e = 0
plot(x1, e0, type = "l", lwd = 2, col = oi["grey"], ylim = c(-0.3, 1.8), xlab = "position",
     ylab = "error", main = "1. smooth  2. coarsen")
lines(x1, e, lwd = 3, col = oi["blue"])
ic <- seq(2, m1, by = 4); points(x1[ic], e[ic], pch = 19, cex = 1.4, col = oi["red"])
legend("topright", bty = "n", cex = 1.15, lwd = c(2, 3, NA), pch = c(NA, NA, 19),
       col = oi[c("grey", "blue", "red")],
       legend = c("error now", "after 5 smoothing sweeps", "smooth error on a coarse grid"))
# right: the V-cycle
par(mar = c(1, 1, 3, 1))
plot(NA, xlim = c(0, 10), ylim = c(0, 5.2), axes = FALSE, xlab = "", ylab = "",
     main = "3. one V-cycle", cex.main = 1.5)
lvl <- c(4.4, 3.2, 2.0, 0.8); xs_ <- c(1, 2.5, 4, 5.5, 7, 8.5)
labs <- c("64 x 64", "32 x 32", "16 x 16", "8 x 8")
xv <- c(1, 2.6, 4.2, 5.8, 7.4, 9.0, 10.6) * 0.92
yv <- c(lvl, rev(lvl)[-1])
lines(xv, yv, lwd = 3, col = oi["grey"])
points(xv, yv, pch = 21, bg = c(rep(oi["blue"], 3), oi["red"], rep(oi["green"], 3)), cex = 3.2)
text(0.3, lvl, labs, cex = 1.15, adj = 1, xpd = NA)
text(2.2, 4.9, "smooth, restrict", cex = 1.2, col = oi["blue"])
text(5.3, 0.2, "solve directly (tiny)", cex = 1.2, col = oi["red"])
text(8.4, 4.9, "interpolate, smooth", cex = 1.2, col = oi["green"])
dev.off()

cat("done\n")
