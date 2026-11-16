# Figures for the least-squares-qr/ unit ("Least squares: normal equations vs QR").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd least-squares-qr && Rscript make-least-squares-qr-figs.R
#
# Base R graphics only; writes SVG files fig01-...svg to fig07-...svg.
# The data are the qPCR and BCA standard curves used on the slides and in the
# exercise (simulated once, then rounded and typed in as fixed numbers).

## ---- the running examples (same numbers as on the slides) -------------------

logN <- rep(1:7, each = 3)                 # log10 copies, 10-fold dilutions, triplicates
cq <- c(34.52, 34.69, 34.48, 31.62, 31.37, 31.14, 28.05, 28.10, 28.07,
        24.54, 24.90, 24.68, 21.13, 20.81, 21.47, 17.89, 17.90, 18.09,
        14.71, 14.67, 14.73)
conc <- rep(c(0, 25, 125, 250, 500, 750, 1000, 1500, 2000), each = 2)  # BSA, ug/mL
a562 <- c(0.087, 0.103, 0.150, 0.109, 0.226, 0.229, 0.358, 0.344, 0.601,
          0.569, 0.780, 0.789, 0.951, 0.942, 1.300, 1.238, 1.545, 1.532)
eps <- .Machine$double.eps

## ---- Gram-Schmidt, classical and modified (same code as on the slides) -------

gs <- function(A, modified = FALSE) {
  n <- ncol(A); Q <- A; R <- matrix(0, n, n)
  for (j in 1:n) {
    v <- A[, j]
    for (i in seq_len(j - 1)) {
      R[i, j] <- sum(Q[, i] * (if (modified) v else A[, j]))
      v <- v - R[i, j] * Q[, i]
    }
    R[j, j] <- sqrt(sum(v^2)); Q[, j] <- v / R[j, j]
  }
  list(Q = Q, R = R)
}

## ---- plotting helpers ---------------------------------------------------------

col_a <- "#0072B2"   # blue
col_b <- "#D55E00"   # vermillion
col_c <- "#009E73"   # bluish green
col_d <- "#CC79A7"   # reddish purple
col_ref <- "grey45"

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

## ---- fig01: qPCR standard curve -------------------------------------------------

fit_q <- lm(cq ~ logN)
b0 <- coef(fit_q)[1]; b1 <- coef(fit_q)[2]
fig("fig01-qpcr-curve.svg", expr = {
  plot(logN, cq, pch = 19, col = col_a, cex = 1.3,
       xlab = expression(log[10] ~ "copies in the standard"), ylab = "Cq (cycles)",
       xlim = c(0.7, 7.3), ylim = c(13, 36))
  abline(fit_q, lwd = 2.5)
  x_unk <- (24 - b0) / b1
  segments(0.5, 24, x_unk, 24, lty = 2, lwd = 2, col = col_b)
  arrows(x_unk, 24, x_unk, 13, lty = 2, lwd = 2, col = col_b, length = 0.12)
  text(0.8, 24.9, "unknown: Cq = 24.0", adj = 0, col = col_b, cex = 1.3)
  text(x_unk + 0.1, 14.5, sprintf("%.2f, about 15 600 copies", x_unk),
       adj = 0, col = col_b, cex = 1.3)
  text(7.2, 34, sprintf("Cq = %.2f - %.3f log10 N\nefficiency %.1f %%", b0, -b1,
                        100 * (10^(-1 / b1) - 1)), adj = 1, cex = 1.3)
})

## ---- fig02: least squares as a projection ---------------------------------------

fig("fig02-projection.svg", w = 8, h = 5.2, expr = {
  par(mar = c(0.5, 0.5, 0.5, 0.5))
  plot(NA, xlim = c(-0.5, 10), ylim = c(-0.5, 6.5), asp = 1, axes = FALSE,
       xlab = "", ylab = "")
  P <- rbind(c(0, 0), c(7, 0), c(9.5, 2.2), c(2.5, 2.2))     # the plane col(A)
  polygon(P, col = adjustcolor(col_a, 0.15), border = col_a, lwd = 1.5)
  text(8.9, 0.6, "all vectors Ax\n(column space of A)", col = col_a, cex = 1.2, adj = 1)
  o <- c(0.8, 0.4); p <- c(5.6, 1.3); b <- c(5.6, 5.8)
  arrows(o[1], o[2], b[1], b[2], lwd = 3, length = 0.15)
  arrows(o[1], o[2], p[1], p[2], lwd = 3, col = col_a, length = 0.15)
  arrows(p[1], p[2], b[1], b[2], lwd = 3, col = col_b, length = 0.15)
  lines(c(p[1], p[1] + 0.35, p[1] + 0.35), c(p[2] + 0.35, p[2] + 0.35, p[2] + 0.08 + 0.0),
        lwd = 1.5)
  text(2.7, 3.6, "b (data)", cex = 1.4)
  text(3.4, 0.55, expression(A * hat(x)), col = col_a, cex = 1.5)
  text(5.8, 3.8, expression(r == b - A * hat(x)), col = col_b, cex = 1.4, adj = 0)
  text(6.1, 5.2, expression(A^T * r == 0), col = col_b, cex = 1.4, adj = 0)
})

## ---- fig03: BCA standard curve, linear vs quadratic --------------------------------

fit1 <- lm(a562 ~ conc); fit2 <- lm(a562 ~ poly(conc, 2))
cc <- seq(0, 2000, length.out = 200)
fig("fig03-bca-fits.svg", w = 11, h = 5, expr = {
  par(mfrow = c(1, 2))
  plot(conc, a562, pch = 19, col = col_a, cex = 1.2,
       xlab = expression("BSA (" * mu * "g/mL)"), ylab = "absorbance 562 nm")
  lines(cc, predict(fit1, data.frame(conc = cc)), lwd = 2.5, lty = 2, col = col_b)
  lines(cc, predict(fit2, data.frame(conc = cc)), lwd = 2.5, col = col_c)
  legend("topleft", c("straight line", "quadratic"), lty = c(2, 1), lwd = 2.5,
         col = c(col_b, col_c), bty = "n", cex = 1.2)
  plot(conc, resid(fit1), pch = 17, col = col_b, cex = 1.3, ylim = c(-0.1, 0.15),
       xlab = expression("BSA (" * mu * "g/mL)"), ylab = "residual")
  points(conc, resid(fit2), pch = 19, col = col_c, cex = 1.2)
  abline(h = 0, col = col_ref)
  legend("topright", c("straight line", "quadratic"), pch = c(17, 19),
         col = c(col_b, col_c), bty = "n", cex = 1.2)
})

## ---- fig04: condition number of the Vandermonde matrix ------------------------------

deg <- 1:6
kraw  <- sapply(deg, function(d) kappa(outer(conc, 0:d, "^"), exact = TRUE))
kne   <- kraw^2
kscl  <- sapply(deg, function(d) kappa(outer(conc / 1000, 0:d, "^"), exact = TRUE))
kpoly <- sapply(deg, function(d) kappa(cbind(1, poly(conc, d)), exact = TRUE))
fig("fig04-vandermonde-kappa.svg", expr = {
  plot(deg, kne, log = "y", type = "b", pch = 15, lwd = 2, col = col_b,
       ylim = c(1, 1e42), yaxt = "n", xlab = "polynomial degree",
       ylab = "condition number")
  logaxis(2, seq(0, 40, by = 10))
  lines(deg, kraw, type = "b", pch = 19, lwd = 2, col = col_a)
  lines(deg, kscl, type = "b", pch = 17, lwd = 2, col = col_c)
  lines(deg, kpoly, type = "b", pch = 18, lwd = 2, col = col_d, cex = 1.4)
  abline(h = 1 / eps, lty = 2, col = col_ref, lwd = 2)
  text(1, 1.3e17, expression(1 / epsilon), col = col_ref, cex = 1.3, adj = 0)
  legend("topleft", c(expression(kappa(A^T * A) ~ ", " * mu * "g/mL"),
                      expression(kappa(A) ~ ", " * mu * "g/mL"),
                      expression(kappa(A) ~ ", mg/mL"),
                      expression(kappa(A) ~ ", poly()")),
         col = c(col_b, col_a, col_c, col_d), pch = c(15, 19, 17, 18),
         lwd = 2, bty = "n", cex = 1.2)
})

## ---- fig05: error against kappa, normal equations vs QR -----------------------------

tt <- seq(0, 1, length.out = 50)
res <- t(sapply(1:12, function(d) {
  A <- outer(tt, 0:d, "^"); xt <- rep(1, d + 1); b <- A %*% xt
  ne <- tryCatch(solve(crossprod(A), crossprod(A, b)), error = function(e) NA)
  q  <- qr.solve(A, b)
  err <- function(x) sqrt(sum((x - xt)^2)) / sqrt(sum(xt^2))
  c(kappa(A, exact = TRUE), err(ne), err(q))
}))
fig("fig05-error-vs-kappa.svg", expr = {
  plot(res[, 1], res[, 2], log = "xy", pch = 15, col = col_b, cex = 1.4,
       xlim = c(1, 1e9), ylim = c(1e-17, 1), xaxt = "n", yaxt = "n",
       xlab = expression(kappa(A)), ylab = "relative error in x")
  logaxis(1, seq(0, 9, by = 1)); logaxis(2, seq(-16, 0, by = 4))
  points(res[, 1], res[, 3], pch = 19, col = col_a, cex = 1.4)
  kk <- 10^seq(0, 9, by = 0.1)
  lines(kk, kk^2 * eps, lty = 2, col = col_b, lwd = 2)
  lines(kk, kk * eps, lty = 2, col = col_a, lwd = 2)
  text(3e4, 3e-4, expression(kappa^2 * epsilon), col = col_b, cex = 1.4)
  text(3e7, 1.5e-10, expression(kappa * epsilon), col = col_a, cex = 1.4, adj = 0)
  legend("topleft", c("normal equations", "QR (qr.solve)"), pch = c(15, 19),
         col = c(col_b, col_a), bty = "n", cex = 1.3)
})

## ---- fig06: a Householder reflection -------------------------------------------------

fig("fig06-householder.svg", w = 7, h = 5.2, expr = {
  par(mar = c(3, 3, 0.5, 0.5))
  plot(NA, xlim = c(-6, 5), ylim = c(-1, 5.5), asp = 1, xlab = "", ylab = "",
       las = 1, cex.axis = 1.2)
  abline(h = 0, v = 0, col = "grey80")
  a <- c(3, 4); v <- a + c(5, 0)                     # v = a + sign(a1) |a| e1
  u <- c(-v[2], v[1]) / sqrt(sum(v^2))               # direction of the mirror
  lines(c(-6, 6) * u[1], c(-6, 6) * u[2], lty = 2, lwd = 2, col = col_ref)
  arrows(0, 0, a[1], a[2], lwd = 3, col = col_a, length = 0.15)
  arrows(0, 0, -5, 0, lwd = 3, col = col_b, length = 0.15)
  arrows(0, 0, v[1] / 2, v[2] / 2, lwd = 2, col = col_c, length = 0.12)
  segments(a[1], a[2], -5, 0, lty = 3, lwd = 1.5)
  text(3.2, 4.3, "a = (3, 4)", col = col_a, cex = 1.4, adj = 0)
  text(-4.6, 0.9, "Ha = (-5, 0)", col = col_b, cex = 1.4)
  text(3.6, 1.2, "direction of\nv = (8, 4)", col = col_c, cex = 1.4)
  text(-2.3, 4.6, "mirror", col = col_ref, cex = 1.3)
})

## ---- fig07: loss of orthogonality ------------------------------------------------------

orth <- function(Q) norm(crossprod(Q) - diag(ncol(Q)), "2")
og <- t(sapply(1:12, function(d) {
  A <- outer(tt, 0:d, "^")
  c(kappa(A, exact = TRUE), orth(gs(A)$Q), orth(gs(A, TRUE)$Q), orth(qr.Q(qr(A))))
}))
fig("fig07-gs-orthogonality.svg", expr = {
  plot(og[, 1], og[, 2], log = "xy", pch = 15, col = col_b, cex = 1.4,
       xlim = c(1, 1e9), ylim = c(1e-17, 10), xaxt = "n", yaxt = "n",
       xlab = expression(kappa(A)), ylab = expression("||" * Q^T * Q - I * "||"))
  logaxis(1, seq(0, 9, by = 1)); logaxis(2, seq(-16, 0, by = 4))
  points(og[, 1], og[, 3], pch = 17, col = col_c, cex = 1.4)
  points(og[, 1], og[, 4], pch = 19, col = col_a, cex = 1.4)
  kk <- 10^seq(0, 9, by = 0.1)
  lines(kk, kk^2 * eps, lty = 2, col = col_b, lwd = 2)
  lines(kk, kk * eps, lty = 2, col = col_c, lwd = 2)
  abline(h = eps, lty = 2, col = col_a, lwd = 2)
  legend("topleft", c("classical Gram-Schmidt", "modified Gram-Schmidt",
                      "Householder (qr)"), pch = c(15, 17, 19),
         col = c(col_b, col_c, col_a), bty = "n", cex = 1.3)
})
