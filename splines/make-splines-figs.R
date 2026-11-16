# Figures for the splines/ unit ("Splines").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd splines && Rscript make-splines-figs.R
#
# Base R graphics (plus pracma for the clamped spline); writes SVG files
# fig01-...svg to fig10-...svg. The data are the noisy OD curve, the gene time
# course and the pseudotime data used on the slides and in the exercise.

library(pracma)   # cubicspline(): clamped end conditions

## ---- the running examples (same numbers as on the slides) -------------------

K <- 1.2; N0 <- 0.05; r <- 0.8                       # logistic OD curve
od   <- function(t) K / (1 + (K / N0 - 1) * exp(-r * t))
rate <- function(t) r * od(t) * (1 - od(t) / K)      # true growth rate N'(t)

t_od <- seq(0, 10, by = 0.25)                        # every 15 min, 41 reads
y_od <- c(0.040, 0.058, 0.076, 0.077, 0.108, 0.127, 0.152, 0.191, 0.200,
          0.263, 0.284, 0.327, 0.382, 0.446, 0.502, 0.556, 0.610, 0.672,
          0.749, 0.794, 0.839, 0.883, 0.934, 0.958, 1.004, 1.032, 1.076,
          1.097, 1.105, 1.111, 1.144, 1.155, 1.163, 1.171, 1.167, 1.182,
          1.193, 1.184, 1.177, 1.197, 1.199)          # od(t) + N(0, 0.01^2), 3 decimals

t_g <- c(0, 0.5, 1, 2, 4, 8, 12, 24)                 # hours after stimulation
e_g <- c(1, 40, 32, 6, 1, 0.5, 0.5, 0.5)             # immediate-early gene, TPM

## ---- plotting helpers ---------------------------------------------------------

col_a <- "#0072B2"   # blue
col_b <- "#D55E00"   # vermillion
col_c <- "#009E73"   # bluish green
col_d <- "#CC79A7"   # reddish purple
col_e <- "#E69F00"   # orange
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
tt <- seq(0, 10, by = 0.01)

## ---- fig01: noisy OD and its finite-difference growth rate ---------------------

n <- length(t_od)
cd <- (y_od[3:n] - y_od[1:(n - 2)]) / 0.5             # central differences
fig("fig01-noisy-od.svg", w = 12, h = 5, expr = {
  par(mfrow = c(1, 2))
  plot(t_od, y_od, pch = 19, col = col_a, xlab = "time (h)", ylab = "OD600")
  lines(tt, od(tt), lwd = 2, col = col_ref, lty = 2)
  legend("topleft", c("plate reader (noise sd 0.01)", "true logistic"),
         pch = c(19, NA), lty = c(NA, 2), lwd = 2, col = c(col_a, col_ref),
         bty = "n", cex = 1.2)
  plot(t_od[2:(n - 1)], cd, pch = 19, col = col_b, ylim = c(-0.05, 0.3),
       xlab = "time (h)", ylab = "growth rate (OD/h)")
  lines(tt, rate(tt), lwd = 2.5, col = col_ref, lty = 2)
  abline(h = 0, col = "grey80")
  legend("topright", c("central difference, h = 0.25 h", "true N'(t)"),
         pch = c(19, NA), lty = c(NA, 2), lwd = 2, col = c(col_b, col_ref),
         bty = "n", cex = 1.2)
})

## ---- fig02: one polynomial vs pieces (gene time course) ------------------------

tg_f <- seq(0, 24, by = 0.01)
V <- outer(t_g, 0:7, "^")
cpoly <- solve(V, e_g)
ppoly <- outer(tg_f, 0:7, "^") %*% cpoly
fig("fig02-local-vs-global.svg", w = 10, h = 5.2, expr = {
  plot(t_g, e_g, pch = 19, cex = 1.5, ylim = c(-15, 50), xlab = "time after stimulation (h)",
       ylab = "expression (TPM)", type = "n")
  abline(h = 0, col = "grey80")
  lines(tg_f, ppoly, col = col_d, lwd = 2)
  lines(t_g, e_g, col = col_c, lwd = 2.5)
  lines(tg_f, splinefun(t_g, e_g, method = "natural")(tg_f), col = col_a, lwd = 2.5)
  points(t_g, e_g, pch = 19, cex = 1.5)
  legend("topright", c("degree-7 polynomial", "piecewise linear", "natural cubic spline"),
         col = c(col_d, col_c, col_a), lwd = 2.5, bty = "n", cex = 1.25)
})

## ---- fig03: error vs spacing ---------------------------------------------------

hs <- 10 / c(5, 10, 20, 40, 80, 160)
err <- t(sapply(hs, function(h) {
  x <- seq(0, 10, by = h); y <- od(x)
  c(lin = max(abs(approx(x, y, tt)$y - od(tt))),
    nat = max(abs(splinefun(x, y, method = "natural")(tt) - od(tt))),
    fmm = max(abs(splinefun(x, y)(tt) - od(tt))),
    der = max(abs(splinefun(x, y)(tt, deriv = 1) - rate(tt))))
}))
fig("fig03-error-orders.svg", w = 9, h = 5.4, expr = {
  matplot(hs, err, log = "xy", type = "b", pch = 19, lty = 1, lwd = 2,
          col = c(col_c, col_b, col_a, col_d), xaxt = "n", yaxt = "n",
          xlab = "spacing h (hours)", ylab = "max error", ylim = c(1e-9, 0.1))
  axis(1, at = c(0.0625, 0.125, 0.25, 0.5, 1, 2),
       labels = c("1/16", "1/8", "1/4", "1/2", "1", "2"))
  logaxis(2, seq(-9, -1, by = 2))
  legend("bottomright", c("piecewise linear (slope 2)", "natural spline (slope 2)",
                          "spline, R default 'fmm' (slope 4)", "its derivative (slope 3)"),
         col = c(col_c, col_b, col_a, col_d), lwd = 2, pch = 19, bty = "n", cex = 1.2)
})

## ---- fig04: boundary conditions ------------------------------------------------

# general cubic spline through (x, y) by a dense solve for the second derivatives
spline_bc <- function(x, y, bc = c("natural", "notaknot", "clamped"), d1 = c(0, 0)) {
  bc <- match.arg(bc)
  n <- length(x); h <- diff(x); s <- diff(y) / h
  A <- matrix(0, n, n); rhs <- numeric(n)
  for (i in 2:(n - 1)) {
    A[i, (i - 1):(i + 1)] <- c(h[i - 1], 2 * (h[i - 1] + h[i]), h[i])
    rhs[i] <- 6 * (s[i] - s[i - 1])
  }
  if (bc == "natural") { A[1, 1] <- 1; A[n, n] <- 1 }
  if (bc == "clamped") {
    A[1, 1:2] <- c(2 * h[1], h[1]);             rhs[1] <- 6 * (s[1] - d1[1])
    A[n, (n - 1):n] <- c(h[n - 1], 2 * h[n - 1]); rhs[n] <- 6 * (d1[2] - s[n - 1])
  }
  if (bc == "notaknot") {     # third derivative continuous at x[2] and x[n-1]
    A[1, 1:3] <- c(h[2], -(h[1] + h[2]), h[1])
    A[n, (n - 2):n] <- c(h[n - 1], -(h[n - 2] + h[n - 1]), h[n - 2])
  }
  M <- solve(A, rhs)
  function(t) {
    i <- findInterval(t, x, all.inside = TRUE)
    hh <- x[i + 1] - x[i]; a <- (x[i + 1] - t) / hh; b <- (t - x[i]) / hh
    a * y[i] + b * y[i + 1] + ((a^3 - a) * M[i] + (b^3 - b) * M[i + 1]) * hh^2 / 6
  }
}
x4 <- 0:10; y4 <- od(x4)
s_nat <- spline_bc(x4, y4, "natural")
s_nak <- spline_bc(x4, y4, "notaknot")
s_cla <- spline_bc(x4, y4, "clamped", d1 = rate(c(0, 10)))
stopifnot(abs(s_cla(0.5) - cubicspline(x4, y4, 0.5, endp2nd = TRUE, der = rate(c(0, 10)))) < 1e-12)
fig("fig04-boundary.svg", w = 10, h = 5.2, expr = {
  plot(tt, 1e3 * (s_nat(tt) - od(tt)), type = "l", lwd = 2.5, col = col_b,
       ylim = c(-1.4, 1.4), xlab = "time (h)", ylab = expression("spline - true OD  (" * 10^-3 * ")"))
  abline(h = 0, col = "grey80"); abline(v = x4, col = "grey90")
  lines(tt, 1e3 * (s_nak(tt) - od(tt)), lwd = 2.5, col = col_a)
  lines(tt, 1e3 * (s_cla(tt) - od(tt)), lwd = 2.5, col = col_c, lty = 2)
  legend("bottomright", c("natural: S'' = 0 at the ends", "not-a-knot", "clamped: S' = true N' at the ends"),
         col = c(col_b, col_a, col_c), lty = c(1, 1, 2), lwd = 2.5, bty = "n", cex = 1.2)
})

## ---- fig05: overshoot ------------------------------------------------------------

fig("fig05-overshoot.svg", w = 10, h = 5.2, expr = {
  plot(t_g, e_g, type = "n", xlim = c(0, 13), ylim = c(-3, 45),
       xlab = "time after stimulation (h)", ylab = "expression (TPM)")
  rect(-1, -5, 30, 0, col = "#F5DCDC", border = NA)
  abline(h = 0, col = "grey60")
  lines(tg_f, splinefun(t_g, e_g)(tg_f), col = col_b, lwd = 2.5)
  lines(tg_f, splinefun(t_g, e_g, method = "monoH.FC")(tg_f), col = col_a, lwd = 2.5)
  lines(tg_f, 2^splinefun(t_g, log2(e_g))(tg_f), col = col_c, lwd = 2.5, lty = 2)
  points(t_g, e_g, pch = 19, cex = 1.5)
  text(3.3, -2.2, "negative expression", col = col_b, cex = 1.2, adj = 0)
  legend("topright", c("cubic spline (R default)", "monotone: method = \"monoH.FC\"",
                       "cubic spline of log2, back-transformed"),
         col = c(col_b, col_a, col_c), lty = c(1, 1, 2), lwd = 2.5, bty = "n", cex = 1.2)
})

## ---- fig06: three smoothing levels -----------------------------------------------

fits <- list(smooth.spline(t_od, y_od, df = 30),
             smooth.spline(t_od, y_od),
             smooth.spline(t_od, y_od, df = 4))
labs <- sprintf(c("too little: df = %.0f", "GCV: df = %.1f", "too much: df = %.0f"),
                sapply(fits, `[[`, "df"))
fig("fig06-lambda.svg", w = 13, h = 6.6, expr = {
  par(mfrow = c(2, 3), mar = c(4.5, 5.6, 2.6, 1))
  for (k in 1:3) {
    plot(t_od, y_od, pch = 19, col = "grey55", xlab = "", ylab = if (k == 1) "OD600" else "",
         main = labs[k], cex.main = 1.6)
    lines(predict(fits[[k]], tt), lwd = 2.5, col = col_a)
  }
  for (k in 1:3) {
    plot(tt, rate(tt), type = "l", lty = 2, lwd = 2.5, col = col_ref, ylim = c(-0.1, 0.35),
         xlab = "time (h)", ylab = if (k == 1) "growth rate (OD/h)" else "")
    lines(predict(fits[[k]], tt, deriv = 1), lwd = 2.5, col = col_b)
  }
})

## ---- fig07: GCV and the true error against df ----------------------------------

lams <- 10^seq(-8, -1.5, by = 0.1)
gc <- t(sapply(lams, function(l) {
  f <- smooth.spline(t_od, y_od, lambda = l)
  d <- predict(f, t_od, deriv = 1)$y
  c(df = f$df, gcv = f$cv.crit, err = sqrt(mean((d - rate(t_od))^2)))
}))
f_gcv <- smooth.spline(t_od, y_od)
fig("fig07-gcv.svg", w = 12, h = 5, expr = {
  par(mfrow = c(1, 2))
  plot(gc[, "df"], gc[, "gcv"], type = "l", lwd = 2.5, col = col_a, log = "xy",
       xlab = "degrees of freedom (df)", ylab = "GCV score", yaxt = "n")
  axis(2, at = c(1e-4, 2e-4, 5e-4, 1e-3), labels = c("1e-4", "2e-4", "5e-4", "1e-3"))
  abline(v = f_gcv$df, lty = 2, col = col_ref)
  text(f_gcv$df * 1.05, max(gc[, "gcv"]) * 0.8, sprintf("min at df = %.1f", f_gcv$df), adj = 0, cex = 1.2)
  plot(gc[, "df"], gc[, "err"], type = "l", lwd = 2.5, col = col_b, log = "xy",
       xlab = "degrees of freedom (df)", ylab = "rms error of growth rate")
  abline(v = f_gcv$df, lty = 2, col = col_ref)
  text(min(gc[, "df"]), max(gc[, "err"]) * 0.6, "  bias\n  dominates", cex = 1.3, adj = 0)
  text(max(gc[, "df"]), max(gc[, "err"]) * 0.6, "noise  \ndominates  ", cex = 1.3, adj = 1)
})

## ---- fig08: growth rate, three ways -------------------------------------------

fig("fig08-growth-rate.svg", w = 10, h = 5.2, expr = {
  plot(t_od[2:(n - 1)], cd, pch = 19, col = "grey60", ylim = c(-0.05, 0.3),
       xlab = "time (h)", ylab = "growth rate (OD/h)")
  abline(h = 0, col = "grey85")
  lines(tt, splinefun(t_od, y_od, method = "natural")(tt, deriv = 1), col = col_d, lwd = 1.5)
  lines(tt, rate(tt), lwd = 3, col = col_ref, lty = 2)
  lines(predict(f_gcv, tt, deriv = 1), lwd = 3, col = col_b)
  legend("topright", c("central differences of the data", "interpolating spline, S'(t)",
                       "smoothing spline (GCV), S'(t)", "true N'(t)"),
         pch = c(19, NA, NA, NA), lty = c(NA, 1, 1, 2), lwd = c(NA, 1.5, 3, 3),
         col = c("grey60", col_d, col_b, col_ref), bty = "n", cex = 1.2)
})

## ---- fig09: pseudotime --------------------------------------------------------

m_pt <- function(t) 1 + 2 * exp(-((t - 0.6) / 0.15)^2)
set.seed(2)
pt3 <- sort(runif(300)); ex3 <- m_pt(pt3) + rnorm(300, sd = 0.5)
f3 <- smooth.spline(pt3, ex3)
g <- seq(0, 1, by = 0.001)
fig("fig09-pseudotime.svg", w = 10, h = 5.2, expr = {
  plot(pt3, ex3, pch = 19, col = adjustcolor(col_a, 0.45), xlab = "pseudotime",
       ylab = "log expression")
  lines(g, m_pt(g), lwd = 2.5, lty = 2, col = col_ref)
  lines(predict(f3, g), lwd = 3, col = col_b)
  legend("topleft", c("300 cells", sprintf("smooth.spline, GCV (df = %.1f)", f3$df), "true mean"),
         pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = c(NA, 3, 2.5),
         col = c(col_a, col_b, col_ref), bty = "n", cex = 1.2)
})

## ---- fig10: correlated noise fools GCV ----------------------------------------

set.seed(3)
y_ar <- od(t_od) + as.numeric(arima.sim(list(ar = 0.9), n = length(t_od), sd = 0.005))
f_ar <- smooth.spline(t_od, y_ar)
fig("fig10-correlated.svg", w = 12, h = 5, expr = {
  par(mfrow = c(1, 2))
  plot(t_od, y_ar, pch = 19, col = "grey55", xlab = "time (h)", ylab = "OD600")
  lines(tt, od(tt), lty = 2, lwd = 2, col = col_ref)
  lines(predict(f_ar, tt), lwd = 2.5, col = col_a)
  legend("topleft", c("drifting reads", sprintf("GCV fit, df = %.1f", f_ar$df), "true"),
         pch = c(19, NA, NA), lty = c(NA, 1, 2), lwd = 2, col = c("grey55", col_a, col_ref),
         bty = "n", cex = 1.2)
  plot(tt, rate(tt), type = "l", lty = 2, lwd = 2.5, col = col_ref, ylim = c(-0.05, 0.35),
       xlab = "time (h)", ylab = "growth rate (OD/h)")
  lines(predict(f_ar, tt, deriv = 1), lwd = 2.5, col = col_b)
  lines(predict(smooth.spline(t_od, y_ar, df = 9.4), tt, deriv = 1), lwd = 2.5, col = col_c)
  legend("topright", c(sprintf("GCV, df = %.1f", f_ar$df), "df = 9.4", "true"),
         lty = c(1, 1, 2), lwd = 2.5, col = c(col_b, col_c, col_ref), bty = "n", cex = 1.2)
})

cat("figures written\n")
