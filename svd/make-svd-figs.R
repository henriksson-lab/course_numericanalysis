# Figures for the svd/ unit ("The SVD, rank deficiency and regularisation").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd svd && Rscript make-svd-figs.R
#
# Base R graphics only (grDevices::svg) plus the nnls package for fig07.
# Colours are Okabe-Ito (colour-blind safe). Makes:
#   fig01-signature.svg       the simulated signature matrix (40 genes x 5 cell types)
#   fig02-rotate-stretch.svg  the SVD of a 2x2 matrix as rotate - stretch - rotate
#   fig03-singular-values.svg singular values of the near and exactly rank-deficient signatures
#   fig04-picard.svg          sigma_i, |u_i^T b| and their ratio (Picard plot)
#   fig05-filter-factors.svg  truncated-SVD and ridge filter factors against sigma
#   fig06-lcurve.svg          L-curve and error against lambda
#   fig07-replicates.svg      estimates over 500 noise replicates: OLS, TSVD, ridge, NNLS

library(nnls)

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

fig <- function(file, w = 9, h = 5.2, expr) {
  svg(file, width = w, height = h)
  par(mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3)
  expr
  dev.off()
}
logaxis <- function(side, at) {
  axis(side, at = 10^at, labels = parse(text = paste0("10^", at)))
}

## ---- the running example (same code as on the slides) ----------------------

g <- 1:40                                  # 40 genes
S <- matrix(10, nrow = 40, ncol = 5,
            dimnames = list(NULL, c("CD4T", "CD8T", "B", "NK", "Mono")))
S[1:16, 1:2] <- 100                        # T-cell genes, high in both T types
S[cbind(17:40, rep(3:5, each = 8))] <- 100 # 8 marker genes each for B, NK, Mono
S <- S * (1 + 0.5 * sin(outer(g, 1:5)))    # gene-to-gene variation
S[, "CD8T"] <- S[, "CD4T"]                 # CD8 = CD4 ...
S[1, "CD8T"] <- 0.98 * S[1, "CD4T"]        # ... except 2% lower on gene 1
S[2, "CD8T"] <- 1.02 * S[2, "CD4T"]        #     and 2% higher on gene 2
p_true <- c(0.30, 0.15, 0.15, 0.10, 0.30)
b0 <- drop(S %*% p_true)
b  <- b0 * (1 + 0.05 * sin(37 * g))        # 5% "noise"
sv <- svd(S)
S0 <- S; S0[, "CD8T"] <- S0[, "CD4T"]      # exactly identical columns

tsvd  <- function(A, b, k) { s <- svd(A)
  drop(s$v[, 1:k, drop = FALSE] %*% (crossprod(s$u[, 1:k, drop = FALSE], b) / s$d[1:k])) }
ridge <- function(A, b, lambda) { s <- svd(A)
  f <- s$d^2 / (s$d^2 + lambda^2)
  drop(s$v %*% (f * crossprod(s$u, b) / s$d)) }

## ---- fig01: signature matrix heatmap ----------------------------------------

fig("fig01-signature.svg", w = 9, h = 5.2, expr = {
  par(mar = c(4.5, 6.5, 1.5, 7))
  pal <- colorRampPalette(c("white", oi["sky"], oi["blue"], "#08306B"))(64)
  z <- log10(t(S))                          # rows of the image = cell types
  image(1:40, 1:5, t(z), col = pal, axes = FALSE, xlab = "gene", ylab = "",
        zlim = range(z))
  axis(1, at = c(1, 8, 16, 24, 32, 40))
  axis(2, at = 1:5, labels = colnames(S), cex.axis = 1.3)
  box()
  abline(v = c(16.5, 24.5, 32.5), col = "grey40", lty = 3)
  # colour key
  par(xpd = NA)
  yk <- seq(1, 5, length.out = 65)
  rect(42.5, yk[-65], 44, yk[-1], col = pal, border = NA)
  rect(42.5, 1, 44, 5)
  lab <- c(5, 10, 30, 100)
  at <- 1 + 4 * (log10(lab) - min(z)) / diff(range(z))
  text(44.4, at, lab, adj = 0, cex = 1.2)
  text(43.3, 5.45, "expression", cex = 1.2)
})

## ---- fig02: rotate - stretch - rotate ----------------------------------------

A <- matrix(c(3, 4, 0, 5), 2, 2)
s2 <- svd(A)
U <- s2$u; V <- s2$v; D <- diag(s2$d)
# choose signs so that v1 and u1 point into the upper right
sg <- sign(V[1, ]); V <- V %*% diag(sg); U <- U %*% diag(sg)
th <- seq(0, 2 * pi, length.out = 200)
circ <- rbind(cos(th), sin(th))
mark <- rbind(c(1, 0), c(0, 1)) # not used; we mark v1, v2 instead
fig("fig02-rotate-stretch.svg", w = 12, h = 4, expr = {
  par(mfrow = c(1, 4), mar = c(1, 1, 3, 1), cex.main = 1.6, font.main = 1)
  panel <- function(P, a1, a2, main, lab1, lab2, lim = 7.2,
                    off1 = c(0, 0.15 * lim), off2 = c(0.15 * lim, 0)) {
    plot(NA, xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1, axes = FALSE,
         xlab = "", ylab = "", main = main)
    abline(h = 0, v = 0, col = "grey85")
    polygon(P[1, ], P[2, ], border = "black", lwd = 2, col = "#56B4E922")
    arrows(0, 0, a1[1], a1[2], col = oi["red"], lwd = 3, length = 0.1)
    arrows(0, 0, a2[1], a2[2], col = oi["blue"], lwd = 3, length = 0.1)
    text(a1[1] + off1[1], a1[2] + off1[2], lab1, col = oi["red"], cex = 1.5)
    text(a2[1] + off2[1], a2[2] + off2[2], lab2, col = oi["blue"], cex = 1.5)
  }
  panel(circ, V[, 1], V[, 2], "x on the unit circle",
        expression(v[1]), expression(v[2]), lim = 2.2)
  P1 <- t(V) %*% circ
  panel(P1, c(1, 0), c(0, 1), expression("rotate: " * V^T * x),
        expression(e[1]), expression(e[2]), lim = 2.2,
        off1 = c(0.2, -0.3), off2 = c(0.3, 0.2))
  P2 <- D %*% P1
  panel(P2, c(s2$d[1], 0), c(0, s2$d[2]), expression("stretch: " * Sigma * V^T * x),
        expression(sigma[1]), expression(sigma[2]),
        off1 = c(-0.3, -1.1), off2 = c(1.1, 0.6))
  P3 <- U %*% P2
  panel(P3, s2$d[1] * U[, 1], s2$d[2] * U[, 2], expression("rotate: " * U * Sigma * V^T * x),
        expression(sigma[1] * u[1]), expression(sigma[2] * u[2]),
        off1 = c(1.6, 0.4), off2 = c(1.6, -0.4))
})

## ---- fig03: singular values, near and exact rank deficiency -----------------

d0 <- svd(S0)$d
fig("fig03-singular-values.svg", expr = {
  plot(1:5, sv$d, log = "y", type = "b", pch = 19, cex = 1.6, lwd = 2,
       col = oi["blue"], ylim = c(1e-15, 2e3), xaxt = "n", yaxt = "n",
       xlab = "index i", ylab = expression(sigma[i]), xlim = c(0.8, 5.4))
  axis(1, at = 1:5); logaxis(2, seq(-15, 3, 3))
  lines(1:5, d0, type = "b", pch = 17, cex = 1.6, lwd = 2, col = oi["red"], lty = 2)
  tol <- 40 * .Machine$double.eps * d0[1]
  abline(h = tol, col = "grey45", lty = 3, lwd = 2)
  text(1, tol * 15, expression("tolerance " * max(m, n) %.% epsilon %.% sigma[1]),
       adj = 0, cex = 1.2, col = "grey30")
  text(4.9, sv$d[5] * 0.15, "2.87", cex = 1.3, col = oi["blue"])
  text(4.9, d0[5] * 0.15, "6.6e-14", cex = 1.3, col = oi["red"])
  legend("bottomleft", c("CD8 differs from CD4 by 2% on 2 genes",
                         "CD8 identical to CD4"),
         col = c(oi["blue"], oi["red"]), pch = c(19, 17), lty = c(1, 2), lwd = 2,
         bty = "n", cex = 1.3)
})

## ---- fig04: Picard plot -----------------------------------------------------

beta  <- abs(drop(crossprod(sv$u, b)))
beta0 <- abs(drop(crossprod(sv$u, b0)))
fig("fig04-picard.svg", expr = {
  plot(1:5, sv$d, log = "y", type = "b", pch = 19, cex = 1.6, lwd = 2,
       col = oi["blue"], ylim = c(5e-3, 2e3), xaxt = "n", yaxt = "n",
       xlab = "index i", ylab = "", xlim = c(0.8, 5.2))
  axis(1, at = 1:5); logaxis(2, -2:3)
  lines(1:5, beta, type = "b", pch = 17, cex = 1.6, lwd = 2, col = oi["red"])
  lines(1:5, beta / sv$d, type = "b", pch = 15, cex = 1.6, lwd = 2, col = oi["green"])
  lines(1:5, beta0 / sv$d, type = "b", pch = 0, cex = 1.6, lwd = 2, lty = 3,
        col = oi["green"])
  legend("bottomleft", c(expression(sigma[i]), expression("|" * u[i]^T * b * "|"),
                         expression("|" * u[i]^T * b * "|" / sigma[i] * "   (noisy b)"),
                         expression("|" * u[i]^T * b[0] * "|" / sigma[i] * "   (noise-free)")),
         col = c(oi["blue"], oi["red"], oi["green"], oi["green"]),
         pch = c(19, 17, 15, 0), lty = c(1, 1, 1, 3), lwd = 2, bty = "n", cex = 1.25)
})

## ---- fig05: filter factors ----------------------------------------------------

sig <- 10^seq(-1, 3.2, length.out = 300)
fig("fig05-filter-factors.svg", expr = {
  plot(sig, sig^2 / (sig^2 + 10^2), log = "x", type = "l", lwd = 3,
       col = oi["blue"], xaxt = "n", ylim = c(0, 1.05),
       xlab = expression("singular value " * sigma), ylab = "filter factor f")
  logaxis(1, -1:3)
  lines(sig, sig^2 / (sig^2 + 1^2), lwd = 3, col = oi["sky"], lty = 2)
  lines(sig, sig^2 / (sig^2 + 100^2), lwd = 3, col = oi["purple"], lty = 4)
  lines(c(0.1, 10, 10, 2000), c(0, 0, 1, 1), lwd = 2.5, col = oi["orange"], lty = 1)
  rug(sv$d, lwd = 3, ticksize = 0.06)
  text(sv$d[5], 0.13, expression(sigma[5]), cex = 1.3)
  text(300, 0.13, expression(sigma[1] * ".." * sigma[4]), cex = 1.3)
  legend("topleft", c(expression("ridge " * lambda == 1), expression("ridge " * lambda == 10),
                      expression("ridge " * lambda == 100), "truncated SVD (cut at 10)"),
         col = c(oi["sky"], oi["blue"], oi["purple"], oi["orange"]),
         lty = c(2, 1, 4, 1), lwd = 3, bty = "n", cex = 1.25)
})

## ---- fig06: L-curve and error against lambda --------------------------------

lam <- 10^seq(-1, 3, by = 0.05)
res <- sol <- err <- numeric(length(lam))
for (i in seq_along(lam)) {
  x <- ridge(S, b, lam[i])
  res[i] <- sqrt(sum((b - S %*% x)^2)); sol[i] <- sqrt(sum(x^2))
  err[i] <- sqrt(sum((x - p_true)^2))
}
fig("fig06-lcurve.svg", w = 12, h = 5.2, expr = {
  par(mfrow = c(1, 2))
  plot(res, sol, log = "xy", type = "l", lwd = 3, col = oi["blue"],
       xlab = "residual norm ||Sx - b||", ylab = "solution norm ||x||",
       ylim = c(0.1, 0.8), xlim = c(5, 220), yaxt = "n")
  axis(2, at = c(0.1, 0.2, 0.4, 0.6, 0.8))
  lab_at <- c(0.1, 3, 10, 100, 1000); lab_pos <- c(4, 4, 2, 1, 2)
  for (i in seq_along(lab_at)) {
    j <- which.min(abs(lam - lab_at[i]))
    points(res[j], sol[j], pch = 19, cex = 1.5, col = oi["red"])
    text(res[j], sol[j], bquote(lambda == .(lab_at[i])), pos = lab_pos[i], cex = 1.25)
  }
  plot(lam, err, log = "x", type = "l", lwd = 3, col = oi["green"], xaxt = "n",
       xlab = expression(lambda), ylab = "error ||x - p_true||", ylim = c(0, 0.45))
  logaxis(1, -1:3)
  abline(h = sqrt(sum((tsvd(S, b, 5) - p_true)^2)), lty = 3, col = "grey45", lwd = 2)
  text(0.12, 0.37, "least squares", adj = 0, cex = 1.2, col = "grey30")
})

## ---- fig07: 500 noise replicates ------------------------------------------------

set.seed(1)
nrep <- 500
est <- array(NA, c(nrep, 4, 3),
             dimnames = list(NULL, c("least\nsquares", "TSVD\nk = 4", "ridge\nλ = 10", "NNLS"),
                             c("CD4", "CD8", "CD4 + CD8")))
for (r in 1:nrep) {
  br <- b0 * (1 + 0.05 * rnorm(40))
  X <- cbind(qr.solve(S, br), tsvd(S, br, 4), ridge(S, br, 10), nnls(S, br)$x)
  est[r, , 1] <- X[1, ]; est[r, , 2] <- X[2, ]; est[r, , 3] <- X[1, ] + X[2, ]
}
fig("fig07-replicates.svg", w = 12, h = 5.2, expr = {
  par(mfrow = c(1, 3), mar = c(5.5, 4.5, 3, 1), cex.main = 1.6, las = 1,
      cex.axis = 1.25, mgp = c(3, 0.8, 0))
  truth <- c(0.30, 0.15, 0.45)
  for (j in 1:3) {
    boxplot(est[, , j], col = c(oi["grey"], oi["orange"], oi["sky"], oi["green"]),
            main = dimnames(est)[[3]][j], ylim = if (j < 3) c(-0.6, 1.05) else c(0.3, 0.6),
            xaxt = "n", outcex = 0.6)
    axis(1, at = 1:4, labels = dimnames(est)[[2]], padj = 0.6, cex.axis = 1.15)
    abline(h = truth[j], col = oi["red"], lwd = 2.5, lty = 2)
    if (j < 3) abline(h = 0, col = "grey60")
  }
})
