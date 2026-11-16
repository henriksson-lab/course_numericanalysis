# Figures for the eigen-algorithms/ unit ("How eigen() works, and big problems").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd eigen-algorithms && Rscript make-eigen-algorithms-figs.R
#
# Base R graphics only (grDevices::svg) plus the Matrix and irlba packages.
# Colours are Okabe-Ito (colour-blind safe). Makes:
#   fig01-pca-cells.svg       PC1 vs PC2 of the simulated 20 000-cell matrix
#   fig02-qr-convergence.svg  subdiagonal entries of the QR algorithm vs (lambda ratio)^k
#   fig03-shifts.svg          last subdiagonal entry: unshifted vs shifted QR
#   fig04-structure.svg       full -> Hessenberg -> triangular; symmetric -> tridiagonal -> diagonal
#   fig05-sensitivity.svg     eigenvalues of a perturbed transit-compartment chain
#   fig06-lanczos.svg         Ritz values against Lanczos step, with and without reorthogonalisation
#   fig07-pca-cost.svg        time and extra memory of full and truncated PCA (measured numbers)
#
# fig07 uses timings measured once (single-threaded OpenBLAS, R 4.6, one desktop
# machine; /usr/bin/time for peak memory). They are typed in below, not re-measured.

suppressMessages({ library(Matrix); library(irlba) })

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999", black = "#000000")

fig <- function(file, w = 9, h = 5.2, expr) {
  svg(file, width = w, height = h)
  par(mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3)
  expr
  dev.off()
}
logaxis <- function(side, pw) axis(side, at = 10^pw, labels = parse(text = paste0("10^", pw)))

## ---- the simulated single-cell matrix (same code as exercise C6) -------------

simulate_cells <- function(n_cells = 20000, n_genes = 2000, n_types = 10, seed = 1) {
  set.seed(seed)
  base <- rgamma(n_genes, shape = 0.3, rate = 2)              # mean count per gene
  fold <- matrix(2^rnorm(n_genes * n_types), n_genes, n_types) # cell-type effects
  type <- sample(n_types, n_cells, replace = TRUE)
  lib  <- rlnorm(n_cells, 0, 0.3)                              # library size factor
  lam  <- base * fold[, type] * rep(lib, each = n_genes)       # genes x cells
  counts <- matrix(rpois(length(lam), lam), n_genes, n_cells)
  list(X = Matrix(t(log1p(counts)), sparse = TRUE), type = type)  # cells x genes
}

## ---- fig01: PCA of the simulated cells ----------------------------------------

sc <- simulate_cells()
X <- sc$X; mu <- colMeans(X)
r <- irlba(X, nv = 10, center = mu)
scores <- r$u %*% diag(r$d)
set.seed(2)
show <- sample(nrow(X), 4000)
tcol <- c(oi[c("orange", "sky", "green", "blue", "red", "purple", "black", "grey")],
          "#8C510A", "#01665E")
fig("fig01-pca-cells.svg", w = 7.5, h = 6, expr = {
  plot(scores[show, 1], scores[show, 2], pch = 16, cex = 0.6,
       col = adjustcolor(tcol[sc$type[show]], 0.7),
       xlab = "PC1 score", ylab = "PC2 score")
  mtext("4000 of 20 000 simulated cells, coloured by cell type", 3, 0.2, cex = 1.2)
})

## ---- fig02: convergence of the unshifted QR algorithm ---------------------------

A <- matrix(c(4, 1, 0, 1, 3, 1, 0, 1, 2), 3)
lam <- eigen(A)$values
K <- 40
sub <- matrix(NA, K + 1, 2)
Ak <- A
sub[1, ] <- abs(c(Ak[2, 1], Ak[3, 2]))
for (k in 1:K) {
  f <- qr(Ak); Ak <- qr.R(f) %*% qr.Q(f)
  sub[k + 1, ] <- abs(c(Ak[2, 1], Ak[3, 2]))
}
fig("fig02-qr-convergence.svg", expr = {
  plot(0:K, sub[, 1], log = "y", type = "o", pch = 19, col = oi["blue"], lwd = 2,
       ylim = c(1e-16, 2), yaxt = "n", xlab = "QR step k", ylab = "|subdiagonal entry|")
  logaxis(2, seq(-16, 0, 4))
  lines(0:K, sub[, 2], type = "o", pch = 17, col = oi["red"], lwd = 2)
  lines(0:K, (lam[2] / lam[1])^(0:K), lty = 2, lwd = 2, col = oi["blue"])
  lines(0:K, (lam[3] / lam[2])^(0:K), lty = 2, lwd = 2, col = oi["red"])
  legend("bottomleft", c(expression("|" * a[21] * "|"), expression("|" * a[32] * "|"),
                         expression((lambda[2] / lambda[1])^k == 0.634^k),
                         expression((lambda[3] / lambda[2])^k == 0.423^k)),
         col = oi[c("blue", "red", "blue", "red")], pch = c(19, 17, NA, NA),
         lty = c(1, 1, 2, 2), lwd = 2, bty = "n", cex = 1.25)
})

## ---- fig03: shifted vs unshifted, last subdiagonal entry --------------------------

last_sub <- function(A, shift, K) {
  out <- numeric(K + 1); n <- nrow(A)
  out[1] <- abs(A[n, n - 1])
  for (k in 1:K) {
    mu <- if (shift) A[n, n] else 0
    f <- qr(A - mu * diag(n)); A <- qr.R(f) %*% qr.Q(f) + mu * diag(n)
    out[k + 1] <- max(abs(A[n, n - 1]), 1e-300)
  }
  out
}
u <- last_sub(A, FALSE, 40)
s <- last_sub(A, TRUE, 4)
fig("fig03-shifts.svg", expr = {
  plot(0:40, u, log = "y", type = "o", pch = 19, col = oi["blue"], lwd = 2,
       ylim = c(1e-17, 2), yaxt = "n", xlab = "QR step k",
       ylab = expression("|" * a[32] * "|"))
  logaxis(2, seq(-16, 0, 4))
  s_plot <- pmax(s, 1e-17)
  lines(0:4, s_plot, type = "o", pch = 17, col = oi["red"], lwd = 2.5)
  text(4.5, 3e-17, "exactly 0", col = oi["red"], adj = 0, cex = 1.2)
  legend("topright", c("no shift: factor 0.42 per step",
                       expression("shift " * mu == a[33] * ": digits triple per step")),
         col = oi[c("blue", "red")], pch = c(19, 17), lwd = 2, bty = "n", cex = 1.25)
})

## ---- fig04: matrix structure along the way ------------------------------------

draw_mat <- function(x0, y0, pat, title) {
  n <- nrow(pat); sz <- 0.32
  for (i in 1:n) for (j in 1:n) {
    col <- if (pat[i, j]) oi["blue"] else "white"
    rect(x0 + (j - 1) * sz, y0 - i * sz, x0 + j * sz, y0 - (i - 1) * sz,
         col = col, border = "grey70")
  }
  text(x0 + n * sz / 2, y0 + 0.3, title, cex = 1.3)
}
n <- 7
full <- matrix(TRUE, n, n)
hess <- outer(1:n, 1:n, function(i, j) i <= j + 1)
tri  <- outer(1:n, 1:n, function(i, j) i <= j)
tridiag <- outer(1:n, 1:n, function(i, j) abs(i - j) <= 1)
dg <- diag(n) == 1
svg("fig04-structure.svg", width = 12, height = 6)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(NA, xlim = c(0, 12.2), ylim = c(-6.2, 0.8), axes = FALSE, xlab = "", ylab = "", asp = 1)
draw_mat(0.2, 0, full, "nonsymmetric A")
draw_mat(4.4, 0, hess, "Hessenberg")
draw_mat(8.6, 0, tri, "triangular (Schur)")
draw_mat(0.2, -3.3, full, "symmetric A")
draw_mat(4.4, -3.3, tridiag, "tridiagonal")
draw_mat(8.6, -3.3, dg, "diagonal")
for (y in c(-1.1, -4.4)) {
  arrows(2.6, y, 4.2, y, lwd = 2.5, length = 0.12)
  arrows(6.8, y, 8.4, y, lwd = 2.5, length = 0.12)
  text(3.4, y - 0.35, "once, O(n^3)", cex = 1.1)
  text(7.6, y - 0.35, "QR steps", cex = 1.1)
}
dev.off()

## ---- fig05: sensitivity of the transit-compartment chain --------------------------

nch <- 10
chain <- diag(-1, nch); chain[cbind(2:nch, 1:(nch - 1))] <- 1
eps_list <- c(1e-14, 1e-10, 1e-6)
ecol <- oi[c("green", "blue", "red")]
fig("fig05-sensitivity.svg", w = 7.5, h = 6.2, expr = {
  par(mar = c(4.8, 5.6, 1.5, 1))
  plot(NA, xlim = c(-1.5, -0.5), ylim = c(-0.42, 0.42), asp = 1,
       xlab = "real part", ylab = "imaginary part")
  abline(h = 0, v = -1, col = "grey80")
  th <- seq(0, 2 * pi, length.out = 200)
  for (i in seq_along(eps_list)) {
    Kp <- chain; Kp[1, nch] <- eps_list[i]
    ev <- eigen(Kp, only.values = TRUE)$values
    rad <- eps_list[i]^(1 / nch)
    lines(-1 + rad * cos(th), rad * sin(th), col = ecol[i], lty = 3)
    points(Re(ev), Im(ev), pch = 19, cex = 1.4, col = ecol[i])
  }
  points(-1, 0, pch = 4, cex = 2, lwd = 3)
  legend("topleft", c(expression(epsilon == 10^-14 * ": moved by 0.04"),
                      expression(epsilon == 10^-10 * ": moved by 0.10"),
                      expression(epsilon == 10^-6 * ":  moved by 0.25"),
                      "unperturbed: all 10 at -1"),
         col = c(ecol, "black"), pch = c(19, 19, 19, 4), bty = "n", cex = 1.1)
})

## ---- fig06: Lanczos Ritz values --------------------------------------------------

lanczos <- function(Afun, n, m, reorth = FALSE) {
  Q <- matrix(0, n, m); alpha <- beta <- numeric(m)
  q <- rep(1, n) / sqrt(n); q_old <- 0; b <- 0
  for (j in 1:m) {
    Q[, j] <- q
    w <- Afun(q) - b * q_old
    alpha[j] <- sum(w * q)
    w <- w - alpha[j] * q
    if (reorth) w <- w - Q[, 1:j] %*% crossprod(Q[, 1:j], w)
    b <- sqrt(sum(w^2)); beta[j] <- b
    q_old <- q; q <- as.vector(w) / b
  }
  list(alpha = alpha, beta = beta)
}
ritz <- function(L, m) {
  T <- diag(L$alpha[1:m], m)
  if (m > 1) T[cbind(2:m, 1:(m - 1))] <- T[cbind(1:(m - 1), 2:m)] <- L$beta[1:(m - 1)]
  eigen(T, symmetric = TRUE, only.values = TRUE)$values
}
spec <- c(10, 9, 8, seq(0, 5, length.out = 197))
M <- 80
Lp <- lanczos(function(v) spec * v, 200, M)
Lr <- lanczos(function(v) spec * v, 200, M, reorth = TRUE)
panel <- function(L, title) {
  plot(NA, xlim = c(1, M), ylim = c(0, 10.6), xlab = "Lanczos step m", ylab = "Ritz values")
  abline(h = c(10, 9, 8), col = "grey80", lwd = 2)
  for (m in 1:M) {
    rv <- ritz(L, m)
    top <- rv[rv > 5.2]
    points(rep(m, length(rv)), rv, pch = 16, cex = 0.35, col = oi["grey"])
    points(rep(m, length(top)), top, pch = 16, cex = 0.6, col = oi["blue"])
  }
  mtext(title, 3, 0.3, cex = 1.25)
}
svg("fig06-lanczos.svg", width = 12, height = 5.2)
par(mfrow = c(1, 2), mar = c(4.8, 5.2, 2.4, 1), las = 1, mgp = c(3.2, 0.8, 0),
    cex.lab = 1.4, cex.axis = 1.3)
panel(Lr, "with reorthogonalisation")
panel(Lp, "plain Lanczos: ghost copies")
dev.off()

## ---- fig07: cost of full vs truncated PCA (measured, see header) -------------------

meth <- c("prcomp\n(dense)", "svd\n(dense)", "eigen\n(covariance)",
          "irlba", "RSpectra\nsvds", "RSpectra\neigs_sym")
secs <- c(21.7, 21.2, 4.7, 0.84, 0.46, 0.42)
mb   <- c(1390, 1125, 1083, 54, 2, 11)   # peak RSS minus 230 MB baseline (R + data)
bcol <- oi[c("red", "red", "orange", "blue", "blue", "blue")]
svg("fig07-pca-cost.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(6.2, 5.6, 2.2, 1), las = 1, mgp = c(3.6, 0.8, 0),
    cex.lab = 1.3, cex.axis = 1.1)
bp <- barplot(secs, col = bcol, ylab = "seconds", names.arg = rep("", 6), ylim = c(0, 25))
text(bp, -1.2, meth, xpd = TRUE, adj = c(0.5, 1), cex = 1.0)
text(bp, secs + 1, format(secs), cex = 1.1)
mtext("time, top 10 PCs", 3, 0.5, cex = 1.3)
bp <- barplot(mb, col = bcol, ylab = "extra memory (MB)", names.arg = rep("", 6), ylim = c(0, 1600))
text(bp, -60, meth, xpd = TRUE, adj = c(0.5, 1), cex = 1.0)
text(bp, mb + 60, format(mb), cex = 1.1)
mtext("peak extra memory", 3, 0.5, cex = 1.3)
dev.off()
