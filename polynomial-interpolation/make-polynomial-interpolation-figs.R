# Figures for the polynomial-interpolation/ unit ("Polynomial interpolation").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd polynomial-interpolation && Rscript make-polynomial-interpolation-figs.R
#
# Base R graphics only; writes SVG files fig01-...svg to fig08-...svg.
# Ground truth is the course's logistic OD curve (K = 1.2, N0 = 0.05,
# r = 0.8 per hour); every interpolant is evaluated with the barycentric
# formula (stable), except the monomial curves in fig08, which show the
# failure of solving the Vandermonde system (forced with solve(..., tol = 0)).

## ---- the running example (same code as on the slides) ----------------------

K <- 1.2; N0 <- 0.05; r <- 0.8          # carrying capacity, start OD, rate per hour
od <- function(t) K / (1 + (K / N0 - 1) * exp(-r * t))

bary_weights <- function(x) {           # scaled to avoid overflow; scaling cancels
  s <- 4 / (max(x) - min(x))
  sapply(seq_along(x), function(j) 1 / prod((x[j] - x[-j]) * s))
}
bary_eval <- function(x, y, w, xx) {
  sapply(xx, function(z) {
    d <- z - x
    if (any(d == 0)) return(y[d == 0][1])
    sum(w * y / d) / sum(w / d)
  })
}
interp <- function(x, y, xx) bary_eval(x, y, bary_weights(x), xx)
cheb <- function(n, a, b) {             # n + 1 Chebyshev points, increasing
  k <- 0:n
  (a + b) / 2 - (b - a) / 2 * cos((2 * k + 1) * pi / (2 * n + 2))
}
nodepoly <- function(xx, nodes) sapply(xx, function(z) prod(z - nodes))

## ---- plotting helpers ---------------------------------------------------------

col_blu <- "#0072B2"   # blue
col_ver <- "#D55E00"   # vermillion
col_grn <- "#009E73"   # bluish green
col_pur <- "#CC79A7"   # reddish purple
col_org <- "#E69F00"   # orange
col_sky <- "#56B4E9"   # sky blue
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

## ---- fig01: five readings, the truth and the degree-4 interpolant -------------

fig("fig01-sparse-od.svg", expr = {
  t <- c(0, 2, 4, 6, 8); y <- od(t)
  tt <- seq(0, 8, length.out = 400)
  plot(tt, od(tt), type = "l", lwd = 3, col = "grey20", ylim = c(0, 1.3),
       xlab = "time (h)", ylab = "OD600")
  lines(tt, interp(t, y, tt), lwd = 3, lty = 2, col = col_blu)
  abline(v = c(3, 5), col = col_ref, lty = 3, lwd = 2)
  points(t, y, pch = 19, cex = 1.8, col = col_ver)
  text(c(3, 5), 0.02, c("RNA 1", "RNA 2"), pos = 4, cex = 1.3, col = col_ref)
  legend("topleft", c("true curve (unknown in practice)", "readings every 2 h",
                      "degree-4 interpolant"),
         col = c("grey20", col_ver, col_blu), lwd = c(3, NA, 3), lty = c(1, NA, 2),
         pch = c(NA, 19, NA), bty = "n", cex = 1.3)
})

## ---- fig02: the five Lagrange basis polynomials -----------------------------

fig("fig02-lagrange-basis.svg", expr = {
  t <- c(0, 2, 4, 6, 8)
  tt <- seq(0, 8, length.out = 400)
  cols <- c(col_blu, col_ver, col_grn, col_pur, col_org)
  plot(NA, xlim = c(0, 8), ylim = c(-0.7, 1.35), xlab = "time (h)",
       ylab = expression(L[i](t)))
  abline(h = c(0, 1), col = "grey80")
  for (i in 1:5) {
    e <- as.numeric(1:5 == i)
    lines(tt, interp(t, e, tt), lwd = 3, col = cols[i])
    points(t, e, pch = 19, cex = 1.4, col = cols[i])
  }
  legend("top", paste0("L[", 0:4, "]  (1 at t = ", t, ")"), col = cols, lwd = 3,
         bty = "n", cex = 1.2, ncol = 3)
})

## ---- fig03: node polynomial, equispaced vs Chebyshev on [0, 48] ---------------

fig("fig03-node-poly.svg", expr = {
  te <- seq(0, 48, by = 4); tc <- cheb(12, 0, 48)
  tt <- seq(0, 48, length.out = 3000)
  we <- abs(nodepoly(tt, te)); wc <- abs(nodepoly(tt, tc))
  plot(tt, we, log = "y", type = "l", lwd = 2.5, col = col_ver, yaxt = "n",
       ylim = c(1e9, 1e16), xlab = "time (h)",
       ylab = expression("|" * w(t) * "|" == prod("|" * t - t[i] * "|")))
  logaxis(2, seq(9, 16, 1))
  lines(tt, wc, lwd = 2.5, col = col_blu)
  rug(te, col = col_ver, lwd = 2, ticksize = 0.04)
  rug(tc, col = col_blu, lwd = 2, ticksize = 0.04, side = 3)
  legend("bottom", c("13 equispaced times (every 4 h)", "13 Chebyshev times"),
         col = c(col_ver, col_blu), lwd = 2.5, bty = "n", cex = 1.3, bg = "white")
})

## ---- fig04: Runge on the 48 h run, equispaced --------------------------------

fig("fig04-runge-od.svg", expr = {
  te <- seq(0, 48, by = 4)
  tt <- seq(0, 48, length.out = 2000)
  pe <- interp(te, od(te), tt)
  plot(tt, od(tt), type = "l", lwd = 3, col = "grey20", ylim = c(-0.25, 1.45),
       xlab = "time (h)", ylab = "OD600")
  abline(h = 0, col = "grey70")
  lines(tt, pe, lwd = 3, col = col_ver)
  points(te, od(te), pch = 19, cex = 1.6, col = col_ver)
  i <- which.min(pe)
  text(tt[i] + 1, pe[i], sprintf("OD = %.2f at %.1f h", pe[i], tt[i]), pos = 4,
       cex = 1.3, col = col_ver)
  legend("bottomright", c("true curve", "degree-12 interpolant, readings every 4 h"),
         col = c("grey20", col_ver), lwd = 3, bty = "n", cex = 1.3)
})

## ---- fig05: same run, Chebyshev times ----------------------------------------

fig("fig05-chebyshev-od.svg", expr = {
  tc <- cheb(12, 0, 48)
  tt <- seq(0, 48, length.out = 2000)
  pc <- interp(tc, od(tc), tt)
  plot(tt, od(tt), type = "l", lwd = 3, col = "grey20", ylim = c(-0.25, 1.45),
       xlab = "time (h)", ylab = "OD600")
  abline(h = 0, col = "grey70")
  lines(tt, pc, lwd = 3, col = col_blu, lty = 2)
  points(tc, od(tc), pch = 19, cex = 1.6, col = col_blu)
  legend("bottomright", c("true curve", "degree-12 interpolant, 13 Chebyshev times"),
         col = c("grey20", col_blu), lwd = 3, lty = c(1, 2), bty = "n", cex = 1.3)
})

## ---- fig06: max error vs degree ----------------------------------------------

maxerr <- function(nodes, a, b) {
  tt <- seq(a, b, length.out = 4001)
  max(abs(interp(nodes, od(nodes), tt) - od(tt)))
}
ns <- 2:40
e8  <- sapply(ns, function(n) maxerr(seq(0, 8, length.out = n + 1), 0, 8))
e48 <- sapply(ns, function(n) maxerr(seq(0, 48, length.out = n + 1), 0, 48))
c48 <- sapply(ns, function(n) maxerr(cheb(n, 0, 48), 0, 48))

fig("fig06-error-vs-n.svg", expr = {
  plot(ns, e48, log = "y", type = "b", pch = 19, col = col_ver, lwd = 2,
       ylim = c(1e-9, 1), yaxt = "n", xlab = "degree n  (n + 1 readings)",
       ylab = "max |error| (OD)")
  logaxis(2, seq(-9, 0, 1))
  lines(ns, c48, type = "b", pch = 17, col = col_blu, lwd = 2)
  lines(ns, e8, type = "b", pch = 15, col = col_grn, lwd = 2)
  legend("bottomleft", c("0-48 h, equispaced", "0-48 h, Chebyshev",
                         "0-8 h, equispaced"),
         col = c(col_ver, col_blu, col_grn), pch = c(19, 17, 15), lwd = 2,
         bty = "n", cex = 1.3)
})

## ---- fig07: interpolation vs least-squares fit on noisy readings -------------

fig("fig07-fit-vs-interp.svg", expr = {
  set.seed(1)
  t12 <- 0:12
  y12 <- od(t12) + rnorm(13, sd = 0.02)
  tt <- seq(0, 12, length.out = 1200)
  f_fit <- lm(y12 ~ poly(t12, 5))
  plot(tt, od(tt), type = "l", lwd = 3, col = "grey20", ylim = c(-0.15, 1.35),
       xlab = "time (h)", ylab = "OD600")
  lines(tt, interp(t12, y12, tt), lwd = 3, col = col_ver)
  lines(tt, predict(f_fit, data.frame(t12 = tt)), lwd = 3, lty = 2, col = col_blu)
  points(t12, y12, pch = 19, cex = 1.5)
  legend("topleft", c("true curve", "noisy readings (sd 0.02)",
                      "interpolant, degree 12", "least squares, degree 5"),
         col = c("grey20", "black", col_ver, col_blu), lwd = c(3, NA, 3, 3),
         lty = c(1, NA, 1, 2), pch = c(NA, 19, NA, NA), bty = "n", cex = 1.3)
})

## ---- fig08: monomial coefficients vs barycentric evaluation ------------------

ns8 <- c(2:20, seq(22, 100, by = 2))
xx <- seq(0, 48, length.out = 2001)
e_bar <- e_mon <- e_raw <- rep(NA, length(ns8))
r_refuse <- NA                                        # first n where solve() refuses
for (j in seq_along(ns8)) {
  n <- ns8[j]; k <- 0:n
  tc <- cheb(n, 0, 48); y <- od(tc)
  w <- (-1)^k * sin((2 * k + 1) * pi / (2 * n + 2))   # Chebyshev weights
  e_bar[j] <- max(abs(bary_eval(tc, y, w, xx) - od(xx)))
  s <- (tc - 24) / 24; ss <- (xx - 24) / 24             # rescaled to [-1, 1]
  cf <- solve(outer(s, 0:n, "^"), y, tol = 0)
  e_mon[j] <- max(abs(outer(ss, 0:n, "^") %*% cf - od(xx)))
  cf2 <- solve(outer(tc, 0:n, "^"), y, tol = 0)       # t in hours, forced
  e_raw[j] <- max(abs(outer(xx, 0:n, "^") %*% cf2 - od(xx)))
  if (is.na(r_refuse) && inherits(try(solve(outer(tc, 0:n, "^"), y), silent = TRUE),
                                  "try-error")) r_refuse <- n
}

fig("fig08-monomial-vs-bary.svg", expr = {
  plot(ns8, e_bar, log = "y", type = "b", pch = 17, col = col_blu, lwd = 2,
       ylim = c(1e-13, 1), yaxt = "n", xlab = "degree n  (Chebyshev times, 0-48 h)",
       ylab = "max |error| (OD)")
  logaxis(2, seq(-13, 0, 2))
  lines(ns8, e_mon, type = "b", pch = 19, col = col_ver, lwd = 2)
  lines(ns8, e_raw, type = "b", pch = 15, col = col_org, lwd = 2, cex = 1.3)
  legend("bottomleft", c("barycentric formula",
                         "solve(Vandermonde), t rescaled to [-1, 1]",
                         "solve(Vandermonde), t in hours"),
         col = c(col_blu, col_ver, col_org), pch = c(17, 19, 15), lwd = 2,
         bty = "n", cex = 1.25)
})

cat("e48 (equispaced 0-48):", signif(e48, 2), "\n")
cat("c48 (Chebyshev 0-48):", signif(c48, 2), "\n")
cat("e8  (equispaced 0-8):", signif(e8, 2), "\n")
cat("solve() on raw Vandermonde refuses from n =", r_refuse, "\n")
cat("raw:", signif(e_raw, 2), "\n")
