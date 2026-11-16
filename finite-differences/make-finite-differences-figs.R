# Figures for the finite-differences/ unit ("Numerical differentiation").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd finite-differences && Rscript make-finite-differences-figs.R
#
# Base R graphics only (plus stats::smooth.spline); writes fig01-...svg to
# fig05-...svg. Everything is computed here, with the same code and data as on
# the slides.

## ---- the running example (same code as on the slides) ----------------------

K <- 1.2; N0 <- 0.05; r <- 0.8          # carrying capacity, start OD, rate per hour
od <- function(t) K / (1 + (K / N0 - 1) * exp(-r * t))
gr <- function(t) r * od(t) * (1 - od(t) / K)   # exact growth rate N'(t)
t0 <- 3
N  <- od(t0); d1 <- gr(t0)
d2 <- r * d1 * (1 - 2 * N / K)                  # N''(3)

# simulated plate reader: every 10 min for 10 h, noise sd 0.005, 3 decimals
set.seed(1)
t  <- seq(0, 10, by = 1/6)
y  <- round(od(t) + rnorm(length(t), sd = 0.005), 3)
dt <- 1/6; n <- length(t)

## ---- plotting helpers ---------------------------------------------------------

col_fwd <- "#0072B2"   # blue
col_cen <- "#D55E00"   # vermillion
col_sec <- "#009E73"   # bluish green
col_ref <- "grey45"
col_pur <- "#CC79A7"   # reddish purple

fig <- function(file, w = 9, h = 5.2, expr, mfrow = c(1, 1)) {
  svg(file, width = w, height = h)
  par(mfrow = mfrow, mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3)
  expr
  dev.off()
}
logaxis <- function(side, at) {
  axis(side, at = 10^at, labels = parse(text = paste0("10^", at)))
}

## ---- fig01: noisy plate-reader data and its naive derivative -----------------

fig("fig01-noisy-od.svg", w = 11, h = 4.8, mfrow = c(1, 2), expr = {
  tt <- seq(0, 10, length.out = 400)
  plot(t, y, pch = 19, cex = 0.9, xlab = "time (h)", ylab = "OD600",
       ylim = c(0, 1.25))
  lines(tt, od(tt), col = col_ref, lwd = 2)
  legend("topleft", c("measured (every 10 min)", "true curve"),
         pch = c(19, NA), lty = c(NA, 1), lwd = c(NA, 2),
         col = c("black", col_ref), bty = "n", cex = 1.2)
  i <- 2:(n - 1)
  cen <- (y[i + 1] - y[i - 1]) / (2 * dt)
  plot(t[i], cen, type = "o", pch = 19, cex = 0.8, col = col_cen,
       xlab = "time (h)", ylab = "growth rate (OD/h)", ylim = c(-0.05, 0.32))
  lines(tt, gr(tt), lwd = 3)
  abline(h = 0, col = "grey80")
  legend("topright", c("central difference", "true N'(t)"),
         col = c(col_cen, "black"), lty = 1, lwd = c(1, 3), pch = c(19, NA),
         bty = "n", cex = 1.2)
})

## ---- fig02: stencils --------------------------------------------------------------

fig("fig02-stencils.svg", w = 9, h = 4.6, expr = {
  par(mar = c(1, 1, 1, 1))
  plot(NA, xlim = c(-3.6, 2.6), ylim = c(0.4, 4.6), axes = FALSE, xlab = "", ylab = "")
  rows <- list(
    list(name = "forward",  s = c(0, 1),     w = c("-1", "1"),       den = "/ h",     col = col_fwd),
    list(name = "backward", s = c(-1, 0),    w = c("-1", "1"),       den = "/ h",     col = col_fwd),
    list(name = "central",  s = c(-1, 1),    w = c("-1", "1"),       den = "/ 2h",    col = col_cen),
    list(name = "second",   s = c(-1, 0, 1), w = c("1", "-2", "1"),  den = "/ h²", col = col_sec))
  for (k in seq_along(rows)) {
    yy <- 5 - k; rr <- rows[[k]]
    segments(-1.5, yy, 1.5, yy, col = "grey70", lwd = 2)
    points(c(-1, 0, 1), rep(yy, 3), pch = 21, bg = "white", col = "grey60", cex = 2.6)
    points(rr$s, rep(yy, length(rr$s)), pch = 19, col = rr$col, cex = 2.6)
    text(rr$s, yy + 0.32, rr$w, cex = 1.4, font = 2)
    text(-3.6, yy, rr$name, adj = 0, cex = 1.45)
    text(1.8, yy, rr$den, adj = 0, cex = 1.45)
  }
  text(c(-1, 0, 1), 0.55, c("t − h", "t", "t + h"), cex = 1.4)
})

## ---- fig03: error against h for three stencils ------------------------------------

eps <- .Machine$double.eps
h   <- 10^-seq(0.5, 12, by = 0.125)
e_f <- abs((od(t0 + h) - od(t0)) / h - d1)
e_c <- abs((od(t0 + h) - od(t0 - h)) / (2 * h) - d1)
e_s <- abs((od(t0 + h) - 2 * od(t0) + od(t0 - h)) / h^2 - d2)
fig("fig03-error-v.svg", w = 9.5, h = 5.2, expr = {
  plot(h, e_f, log = "xy", type = "l", lwd = 3, col = col_fwd, xaxt = "n", yaxt = "n",
       ylim = c(1e-15, 1e4), xlab = "step h", ylab = "absolute error")
  lines(h, e_c, lwd = 3, col = col_cen)
  lines(h, e_s, lwd = 3, col = col_sec)
  logaxis(1, seq(-12, 0, by = 2)); logaxis(2, seq(-15, 3, by = 3))
  abline(v = c(eps^(1/2), eps^(1/3), eps^(1/4)), lty = 3, lwd = 2,
         col = c(col_fwd, col_cen, col_sec))
  legend("topright", c("forward  N'  (slope 1)", "central  N'  (slope 2)",
                       "second-derivative  N''  (slope 2)"),
         col = c(col_fwd, col_cen, col_sec), lwd = 3, bty = "n", cex = 1.15,
         bg = "white")
  text(eps^(1/2), 1e-14, expression(epsilon^{1/2}), pos = 4, cex = 1.3, col = col_fwd)
  text(eps^(1/3), 1e-14, expression(epsilon^{1/3}), pos = 4, cex = 1.3, col = col_cen)
  text(eps^(1/4), 1e-14, expression(epsilon^{1/4}), pos = 4, cex = 1.3, col = col_sec)
})

## ---- fig04: noisy data, error against half-width ---------------------------------

k   <- 1:12
rms <- sapply(k, function(k) {
  j <- (k + 1):(n - k)
  sqrt(mean(((y[j + k] - y[j - k]) / (2 * k * dt) - gr(t[j]))^2))
})
hw <- k * dt
fig("fig04-noise-v.svg", w = 9, h = 5.0, expr = {
  plot(hw, rms, log = "xy", type = "o", pch = 19, lwd = 3, cex = 1.3, col = col_cen,
       xlab = "half-width of the stencil h (hours)", ylab = "rms error of N' (OD/h)",
       ylim = c(0.002, 0.05))
  trunc <- sapply(k, function(k) {           # same, but with exact OD values
    j <- (k + 1):(n - k)
    sqrt(mean(((od(t[j + k]) - od(t[j - k])) / (2 * k * dt) - gr(t[j]))^2))
  })
  lines(hw, 0.005 / (sqrt(2) * hw), lty = 2, lwd = 2, col = col_ref)
  lines(hw, trunc, lty = 3, lwd = 2, col = col_ref)
  text(0.3, 0.012, expression("noise " * sigma / (sqrt(2) * h)), cex = 1.3, col = col_ref, pos = 4)
  text(0.62, 0.0021, "truncation", cex = 1.3, col = col_ref, pos = 4)
  text(dt * 1.08, rms[1] * 1.4, "neighbours (10 min)", adj = 0, cex = 1.15)
})

## ---- fig05: smoothing before differentiating ---------------------------------------

fig("fig05-smooth.svg", w = 9, h = 5.0, expr = {
  tt <- seq(0, 10, length.out = 400)
  i <- 2:(n - 1)
  cen <- (y[i + 1] - y[i - 1]) / (2 * dt)
  k <- 4; j <- (k + 1):(n - k)
  wide <- (y[j + k] - y[j - k]) / (2 * k * dt)
  fit <- smooth.spline(t, y)
  plot(t[i], cen, pch = 19, cex = 0.8, col = "grey65", xlab = "time (h)",
       ylab = "growth rate (OD/h)", ylim = c(-0.05, 0.32))
  lines(tt, gr(tt), lwd = 5)
  lines(t[j], wide, type = "o", pch = 15, cex = 0.8, lwd = 2, col = col_fwd)
  lines(tt, predict(fit, tt, deriv = 1)$y, lwd = 3, col = col_cen, lty = 1)
  abline(h = 0, col = "grey80")
  legend("topright", c("true N'(t)", "central, neighbours", "central, h = 40 min",
                       "smoothing spline"),
         col = c("black", "grey65", col_fwd, col_cen), lwd = c(5, NA, 2, 3),
         pch = c(NA, 19, 15, NA), bty = "n", cex = 1.15)
})
