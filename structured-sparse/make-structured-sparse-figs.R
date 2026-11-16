# Figures for the structured-sparse/ unit ("Structure: symmetric, banded, sparse").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd structured-sparse && Rscript make-structured-sparse-figs.R
#
# Base R graphics plus the Matrix package (ships with R). Colours are
# Okabe-Ito (colour-blind safe). Makes:
#   fig01-structures.svg   three structures: SPD (dense), tridiagonal, sparse
#   fig02-morphogen.svg    steady-state morphogen gradient in a row of cells
#   fig03-timing.svg       dense solve vs Thomas vs sparse solve, time vs n
#   fig04-fill-in.svg      arrow matrix: Cholesky factor with hub first / last
#   fig05-iterative.svg    residual vs iteration for Jacobi and CG
#
# Everything is computed here except fig03, whose timings were measured once
# (R 4.6, a Linux server, OPENBLAS_NUM_THREADS=1). Set MEASURE <- TRUE to
# re-measure on your own machine (about a minute); your numbers will differ,
# the slopes should not.

MEASURE <- FALSE

suppressMessages(library(Matrix))

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- the running example (same code as on the slides) ----------------------

thomas <- function(a, b, c, d) {
  n <- length(d)
  for (i in 2:n) {                 # forward sweep: eliminate a[i]
    w <- a[i] / b[i - 1]
    b[i] <- b[i] - w * c[i - 1]
    d[i] <- d[i] - w * d[i - 1]
  }
  x <- numeric(n)
  x[n] <- d[n] / b[n]
  for (i in (n - 1):1)             # back substitution
    x[i] <- (d[i] - c[i] * x[i + 1]) / b[i]
  x
}

# morphogen in a row of n cells: exchange D with each neighbour, decay k,
# production 1 in cell 1, no flux out of either end
morphogen <- function(n, D = 1, k = 0.01) {
  list(a = c(0, rep(-D, n - 1)), cc = c(rep(-D, n - 1), 0),
       b = c(D + k, rep(2 * D + k, n - 2), D + k), d = c(1, rep(0, n - 1)))
}
as_sparse <- function(m) {
  n <- length(m$d)
  bandSparse(n, k = -1:1, diagonals = list(m$a[-1], m$b, m$cc[-n]))
}

# draw the non-zero pattern of a matrix (row 1 at the top)
spy <- function(M, main, col = "black", zero_col = "white", value = FALSE) {
  M <- as.matrix(M)
  z <- if (value) abs(M) else (abs(M) > 1e-14) * 1
  z <- t(z[nrow(z):1, , drop = FALSE])
  if (value) {
    image(z, col = gray.colors(30, start = 1, end = 0.1), axes = FALSE, useRaster = FALSE)
  } else {
    image(z, col = c(zero_col, col), zlim = c(0, 1), axes = FALSE, useRaster = FALSE)
  }
  box(col = oi["grey"])
  title(main, cex.main = 1.9, font.main = 1, line = 0.8)
}

## ---- fig01: three structures ------------------------------------------------

svg("fig01-structures.svg", width = 12, height = 4.6)
par(mfrow = c(1, 3), mar = c(2.6, 1, 3.2, 1))
n <- 30
spy(0.9^abs(outer(1:n, 1:n, "-")), "symmetric positive definite", value = TRUE)
mtext("covariance, kinship: dense", side = 1, line = 1, cex = 1.3)
spy(as_sparse(morphogen(n)), "tridiagonal", col = oi["blue"])
mtext("row of cells: 3n non-zeros", side = 1, line = 1, cex = 1.3)
set.seed(1)
cnt <- matrix(rpois(30 * 45, 0.1), nrow = 30)
spy(cnt, "sparse", col = oi["red"])
mtext(sprintf("genes x cells: %d %% non-zero", round(100 * mean(cnt > 0))),
      side = 1, line = 1, cex = 1.3)
dev.off()

## ---- fig02: morphogen gradient ----------------------------------------------

n <- 100
m <- morphogen(n)
x <- thomas(m$a, m$b, m$cc, m$d)
svg("fig02-morphogen.svg", width = 9, height = 5.2)
par(mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.6, 0.8, 0),
    cex.lab = 1.5, cex.axis = 1.3)
plot(1:n, x, log = "y", pch = 16, cex = 0.8, col = oi["blue"],
     xlab = "cell number", ylab = "concentration (steady state)",
     ylim = c(5e-4, 20), yaxt = "n")
axis(2, at = 10^(-3:1), labels = c("0.001", "0.01", "0.1", "1", "10"))
lines(1:n, x[1] * exp(-(1:n - 1) / 10), col = oi["orange"], lwd = 2, lty = 2)
legend("topright", c("Thomas solve, 100 cells", "x[1] exp(-(i - 1) / 10)"),
       col = c(oi["blue"], oi["orange"]), pch = c(16, NA), lty = c(NA, 2),
       lwd = c(NA, 2), cex = 1.3, bty = "n")
text(55, 2, "decay length sqrt(D / k) = 10 cells", cex = 1.3)
dev.off()

## ---- fig03: timing ----------------------------------------------------------

nd <- c(500, 1000, 2000, 4000, 8000)
ns <- c(1e4, 3e4, 1e5, 3e5, 1e6, 3e6)
if (MEASURE) {
  tm <- function(f, reps = 3) median(replicate(reps, system.time(f())["elapsed"]))
  dense <- sapply(nd, function(n) {
    m <- morphogen(n); A <- as.matrix(as_sparse(m)); tm(function() solve(A, m$d)) })
  th <- sapply(ns, function(n) {
    m <- morphogen(n); tm(function() thomas(m$a, m$b, m$cc, m$d)) })
  sp <- sapply(ns, function(n) {
    m <- morphogen(n); A <- as_sparse(m); tm(function() solve(A, m$d)) })
  print(list(dense = dense, thomas = th, sparse = sp))
} else {                                   # measured, see header
  dense <- c(0.005, 0.029, 0.144, 0.919, 6.308)
  th    <- c(0.004, 0.015, 0.052, 0.160, 0.551, 1.692)
  sp    <- c(0.001, 0.004, 0.012, 0.040, 0.130, 0.482)
}
svg("fig03-timing.svg", width = 9, height = 5.2)
par(mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.6, 0.8, 0),
    cex.lab = 1.5, cex.axis = 1.3)
plot(NA, log = "xy", xlim = c(400, 4e6), ylim = c(1e-4, 20), xaxt = "n", yaxt = "n",
     xlab = "number of cells n", ylab = "time (seconds)")
axis(1, at = 10^(3:6), labels = c(expression(10^3), expression(10^4),
                                   expression(10^5), expression(10^6)))
axis(2, at = 10^(-4:1), labels = c("0.0001", "0.001", "0.01", "0.1", "1", "10"))
lines(nd, dense, type = "b", pch = 16, lwd = 2, col = oi["red"])
lines(ns, th, type = "b", pch = 17, lwd = 2, col = oi["blue"])
lines(ns, sp, type = "b", pch = 15, lwd = 2, col = oi["green"])
s1 <- coef(lm(log10(dense) ~ log10(nd)))[2]
s2 <- coef(lm(log10(th) ~ log10(ns)))[2]
s3 <- coef(lm(log10(sp[-1]) ~ log10(ns[-1])))[2]
legend("bottomright", bty = "n", cex = 1.25, lwd = 2, pch = c(16, 17, 15),
       col = c(oi["red"], oi["blue"], oi["green"]),
       legend = c(sprintf("dense solve(): slope %.1f", s1),
                  sprintf("Thomas, R loop: slope %.1f", s2),
                  sprintf("sparse solve() (Matrix): slope %.1f", s3)))
dev.off()

## ---- fig04: fill-in for an arrow matrix -------------------------------------

arrow <- function(n, hub) {
  A <- diag(2, n); A[hub, ] <- 1; A[, hub] <- 1; A[hub, hub] <- n
  A
}
n <- 12
A1 <- arrow(n, 1); A2 <- arrow(n, n)
svg("fig04-fill-in.svg", width = 12, height = 4.4)
par(mfrow = c(1, 3), mar = c(1, 1, 3.2, 1))
spy(A1, "A, hub gene first", col = oi["blue"])
spy(t(chol(A1)), "L, hub first: all filled", col = oi["red"])
spy(t(chol(A2)), "L, hub last: no fill-in", col = oi["green"])
dev.off()

## ---- fig05: Jacobi vs conjugate gradient ------------------------------------

n <- 1000
A <- as_sparse(morphogen(n)); bvec <- morphogen(n)$d
jacobi <- function(A, b, iters) {
  dA <- diag(A); x <- numeric(length(b)); res <- numeric(iters)
  for (it in 1:iters) {
    x <- x + (b - A %*% x)[, 1] / dA
    res[it] <- sqrt(sum((b - A %*% x)^2))
  }
  res
}
cg <- function(A, b, iters) {
  x <- numeric(length(b)); r <- b; p <- r; rr <- sum(r^2); res <- numeric(iters)
  for (it in 1:iters) {
    Ap <- (A %*% p)[, 1]
    alpha <- rr / sum(p * Ap)
    x <- x + alpha * p; r <- r - alpha * Ap
    rr_new <- sum(r^2); res[it] <- sqrt(rr_new)
    p <- r + (rr_new / rr) * p; rr <- rr_new
  }
  res
}
rj <- jacobi(A, bvec, 5000); rc <- cg(A, bvec, 400)
cat("Jacobi below 1e-10 at", which(rj < 1e-10)[1],
    "; CG below 1e-10 at", which(rc < 1e-10)[1], "\n")
svg("fig05-iterative.svg", width = 9, height = 5.2)
par(mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.6, 0.8, 0),
    cex.lab = 1.5, cex.axis = 1.3)
plot(NA, log = "y", xlim = c(0, 5000), ylim = c(1e-13, 2), yaxt = "n",
     xlab = "iteration", ylab = "residual  ||b - A x||")
axis(2, at = 10^seq(-12, 0, by = 3),
     labels = c(expression(10^-12), expression(10^-9), expression(10^-6),
                expression(10^-3), "1"))
lines(seq_along(rj), rj, lwd = 2.5, col = oi["blue"])
cc_end <- which(rc < 1e-13)[1]; if (is.na(cc_end)) cc_end <- length(rc)
lines(seq_len(cc_end), rc[seq_len(cc_end)], lwd = 2.5, col = oi["red"])
abline(h = 1e-10, lty = 3, col = oi["grey"])
legend("topright", bty = "n", cex = 1.3, lwd = 2.5, col = c(oi["blue"], oi["red"]),
       legend = c("Jacobi", "conjugate gradient"))
dev.off()
