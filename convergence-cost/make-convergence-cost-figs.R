# Figures for the convergence-cost/ unit ("Convergence and cost").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd convergence-cost && Rscript make-convergence-cost-figs.R
#
# Base R graphics only; writes SVG files fig01-...svg to fig07-...svg.
# Everything is computed here except fig06, whose timings were measured once
# with the code on the "A doubling experiment" slide (R 4.6, a Linux server).
# Set MEASURE <- TRUE to re-measure on your own machine (takes ~20 s); your
# numbers will differ, the slope should not.

MEASURE <- FALSE

## ---- the running example (same code as on the slides) ----------------------

K <- 1.2; N0 <- 0.05; r <- 0.8          # carrying capacity, start OD, rate per hour
od <- function(t) K / (1 + (K / N0 - 1) * exp(-r * t))
t0 <- 3
exact <- r * od(t0) * (1 - od(t0) / K)  # logistic: N' = r N (1 - N/K)
f2 <- r * exact * (1 - 2 * od(t0) / K)  # second derivative at t0
eps <- .Machine$double.eps

h   <- 10^-seq(1, 15, by = 0.25)
fwd <- abs((od(t0 + h) - od(t0)) / h - exact)
cen <- abs((od(t0 + h) - od(t0 - h)) / (2 * h) - exact)

## ---- plotting helpers ---------------------------------------------------------

col_fwd <- "#0072B2"   # blue
col_cen <- "#D55E00"   # vermillion
col_ref <- "grey45"
col_grn <- "#009E73"   # bluish green
col_pur <- "#CC79A7"   # reddish purple

fig <- function(file, w = 9, h = 5.2, expr) {
  svg(file, width = w, height = h)
  par(mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3)
  expr
  dev.off()
}

# log10 axis with labels 10^k
logaxis <- function(side, at) {
  axis(side, at = 10^at, labels = parse(text = paste0("10^", at)))
}

## ---- fig01: OD curve, tangent and secant ---------------------------------------

fig("fig01-od-curve.svg", expr = {
  tt <- seq(0, 10, length.out = 400)
  plot(tt, od(tt), type = "l", lwd = 3, xlab = "time (h)", ylab = "OD600",
       ylim = c(0, 1.25))
  hh <- 1.5
  slope_sec <- (od(t0 + hh) - od(t0)) / hh
  ts <- c(0.5, 7)
  lines(ts, od(t0) + exact * (ts - t0), col = col_grn, lwd = 3)
  lines(ts, od(t0) + slope_sec * (ts - t0), col = col_cen, lwd = 3, lty = 2)
  points(c(t0, t0 + hh), od(c(t0, t0 + hh)), pch = 19, cex = 1.6)
  segments(t0, od(t0), t0 + hh, od(t0), col = "grey40")
  segments(t0 + hh, od(t0), t0 + hh, od(t0 + hh), col = "grey40")
  text(t0 + hh / 2, od(t0) - 0.05, "h", cex = 1.4)
  text(t0 + hh + 0.1, (od(t0) + od(t0 + hh)) / 2, "f(t+h) - f(t)", adj = 0, cex = 1.3)
  legend("bottomright", c("true slope f'(t) = growth rate",
                          "forward difference (f(t+h) - f(t)) / h"),
         col = c(col_grn, col_cen), lwd = 3, lty = c(1, 2), bty = "n", cex = 1.3)
})

## ---- fig02: the V-shaped error of the forward difference --------------------

fig("fig02-error-v.svg", expr = {
  plot(h, fwd, log = "xy", type = "b", pch = 19, col = col_fwd, lwd = 2,
       xaxt = "n", yaxt = "n", ylim = c(1e-10, 1),
       xlab = "step h", ylab = "|error|")
  logaxis(1, seq(-15, -1, 2)); logaxis(2, seq(-10, 0, 2))
  lines(h, abs(f2) / 2 * h, lty = 2, lwd = 2, col = col_ref)
  lines(h, eps * od(t0) / h, lty = 3, lwd = 2.5, col = col_ref)
  hopt <- sqrt(2 * eps * od(t0) / abs(f2))
  abline(v = hopt, col = col_grn, lwd = 2)
  text(1e-4, 1e-3, "truncation\n~ h", cex = 1.3, col = col_ref)
  text(1e-13, 3e-7, "rounding\n~ eps / h", cex = 1.3, col = col_ref)
  text(hopt * 1.6, 1e-1, expression(h %~~% sqrt(epsilon)), adj = 0, cex = 1.4, col = col_grn)
})

## ---- fig03: forward vs central difference, slopes 1 and 2 -------------------

fig("fig03-orders.svg", expr = {
  plot(h, fwd, log = "xy", type = "b", pch = 19, col = col_fwd, lwd = 2,
       xaxt = "n", yaxt = "n", ylim = c(1e-13, 1),
       xlab = "step h", ylab = "|error|")
  logaxis(1, seq(-15, -1, 2)); logaxis(2, seq(-12, 0, 2))
  lines(h, cen, type = "b", pch = 17, col = col_cen, lwd = 2)
  # slope triangles in the truncation region, drawn beside each curve
  tri <- function(x1, x2, y1, p, col, upper) {
    y2 <- y1 * (x2 / x1)^p
    lines(c(x1, x2), c(y1, y2), col = col, lwd = 2)
    if (upper) {
      lines(c(x1, x1, x2), c(y1, y2, y2), col = col, lwd = 1.5, lty = 2)
      text(x1 / 1.4, sqrt(y1 * y2), paste("slope", p), adj = 1, col = col, cex = 1.3)
    } else {
      lines(c(x1, x2, x2), c(y1, y1, y2), col = col, lwd = 1.5, lty = 2)
      text(x2 * 1.4, sqrt(y1 * y2), paste("slope", p), adj = 0, col = col, cex = 1.3)
    }
  }
  tri(1e-6, 1e-4, 3e-6, 1, col_fwd, upper = TRUE)
  tri(1e-4, 1e-2, 7e-13, 2, col_cen, upper = FALSE)
  legend("bottomleft", c("forward  (f(t+h) - f(t)) / h",
                         "central  (f(t+h) - f(t-h)) / 2h"),
         col = c(col_fwd, col_cen), pch = c(19, 17), lwd = 2, bty = "n", cex = 1.3)
})

## ---- fig04: linear vs quadratic convergence (sqrt 2) --------------------------

fig("fig04-iterations.svg", expr = {
  nstep <- 12
  lo <- 1; hi <- 2; x <- 1
  e_bis <- e_new <- numeric(nstep)
  for (k in 1:nstep) {
    mid <- (lo + hi) / 2
    if (mid^2 > 2) hi <- mid else lo <- mid
    e_bis[k] <- abs(mid - sqrt(2))
    x <- (x + 2 / x) / 2
    e_new[k] <- abs(x - sqrt(2))
  }
  e_new[e_new == 0] <- NA
  plot(1:nstep, e_bis, log = "y", type = "b", pch = 19, col = col_fwd, lwd = 2,
       ylim = c(1e-16, 1), yaxt = "n", xlab = "step k", ylab = "|error|")
  logaxis(2, seq(-16, 0, 4))
  lines(1:nstep, 2^-(1:nstep), lty = 2, col = col_fwd, lwd = 1.5)
  lines(1:nstep, e_new, type = "b", pch = 17, col = col_cen, lwd = 2)
  abline(h = eps, col = col_ref, lty = 3, lwd = 2)
  text(12, eps * 8, "machine epsilon", adj = 1, col = col_ref, cex = 1.2)
  legend(5.3, 1e-5, c("halving (bisection): linear",
                         "bisection bound 2^-k",
                         "Newton: quadratic"),
         col = c(col_fwd, col_fwd, col_cen), pch = c(19, NA, 17),
         lty = c(1, 2, 1), lwd = 2, bty = "n", cex = 1.3)
})

## ---- fig05: operation counts vs n ----------------------------------------------

fig("fig05-op-counts.svg", expr = {
  n <- 10^seq(1, 5, by = 0.05)
  flops <- 1e9                                # assumed: 1e9 operations per second
  plot(n, 2 * n^3 / flops, log = "xy", type = "l", lwd = 3, col = col_cen,
       ylim = c(1e-8, 1e7), xaxt = "n", yaxt = "n",
       xlab = "n", ylab = "seconds at 1e9 ops / s")
  logaxis(1, 1:5); logaxis(2, seq(-8, 6, 2))
  lines(n, 2 * n^2 / flops, lwd = 3, col = col_pur)
  lines(n, 2 * n / flops, lwd = 3, col = col_fwd)
  abline(h = c(1, 3600, 86400), col = col_ref, lty = 3)
  text(3000, c(1, 3600, 86400) * 2.2, c("1 second", "1 hour", "1 day"),
       adj = 0, col = col_ref, cex = 1.2)
  abline(v = c(100, 20000), col = col_grn, lwd = 2)
  text(c(100, 20000) * 1.15, 1e-7, c("100", "20 000"), adj = 0, col = col_grn, cex = 1.3)
  legend("topleft", inset = c(0.08, 0), c("matrix-matrix  2n^3", "matrix-vector  2n^2", "dot product  2n"),
         col = c(col_cen, col_pur, col_fwd), lwd = 3, bty = "n", cex = 1.2)
})

## ---- fig06: timing of cor() vs number of genes ----------------------------------

genes <- c(500, 1000, 2000, 4000)
if (MEASURE) {
  set.seed(1)
  cells <- 1000
  secs  <- numeric(length(genes))
  for (i in seq_along(genes)) {
    X <- matrix(rpois(cells * genes[i], 2), nrow = cells)   # cells x genes
    secs[i] <- system.time(cor(X))["elapsed"]
  }
  print(secs)
} else {
  secs <- c(0.167, 0.654, 2.594, 10.946)   # measured, see header
}
fit <- lm(log10(secs) ~ log10(genes))

fig("fig06-cor-timing.svg", expr = {
  gx <- c(300, 30000)
  plot(genes, secs, log = "xy", pch = 19, cex = 1.6, col = col_fwd,
       xlim = gx, ylim = c(0.05, 1000), xaxt = "n", yaxt = "n",
       xlab = "number of genes (1000 cells)", ylab = "seconds for cor()")
  axis(1, at = c(300, 1000, 3000, 10000, 20000),
       labels = c("300", "1000", "3000", "10 000", "20 000"))
  axis(2, at = c(0.1, 1, 10, 100, 1000), labels = c("0.1", "1", "10", "100", "1000"))
  gg <- 10^seq(log10(gx[1]), log10(gx[2]), length.out = 50)
  lines(gg, 10^predict(fit, data.frame(genes = gg)), lty = 2, lwd = 2, col = col_ref)
  pred <- 10^predict(fit, data.frame(genes = 20000))
  points(20000, pred, pch = 1, cex = 2, lwd = 2, col = col_cen)
  text(20000, pred * 2.5, sprintf("20 000 genes:\n~%.0f s predicted", pred),
       col = col_cen, cex = 1.2, adj = 1)
  text(400, 200, sprintf("fitted slope = %.2f", coef(fit)[2]), adj = 0, cex = 1.4)
})

## ---- fig07: memory -----------------------------------------------------------

fig("fig07-memory.svg", w = 9, h = 5, expr = {
  par(mar = c(4.8, 17, 1.5, 1))
  gb <- c("dense, 100 x 100 000" = 100 * 1e5 * 8 / 1e9,
          "correlations, 20 000 x 20 000" = 20000^2 * 8 / 1e9,
          "sparse 10 %, 20 000 x 100 000" = 20000 * 1e5 * 0.1 * 12 / 1e9,
          "dense, 20 000 x 100 000" = 20000 * 1e5 * 8 / 1e9)
  b <- barplot(gb, horiz = TRUE, names.arg = "", xlim = c(0, 19),
               col = c(col_grn, col_pur, col_fwd, col_cen), border = NA,
               xlab = "memory (GB)")
  axis(2, at = b, labels = names(gb), tick = FALSE, cex.axis = 1.3)
  text(gb + 0.2, b, sprintf("%.2g GB", gb), adj = 0, cex = 1.2)
  abline(v = 16, lty = 2, lwd = 2, col = col_ref)
  mtext("16 GB laptop", side = 3, at = 16, line = 0.1, col = col_ref, cex = 1.1)
})

## ---- numbers used in the text --------------------------------------------------

cat("od(t0) =", od(t0), " exact slope =", exact, " f'' =", f2, "\n")
cat("best h (forward) =", h[which.min(fwd)], " error", min(fwd), "\n")
cat("hopt formula =", sqrt(2 * eps * od(t0) / abs(f2)), "\n")
cat("best h (central) =", h[which.min(cen)], " error", min(cen), "\n")
cat("cor() slope =", coef(fit)[2], " pred 20000 genes =",
    10^predict(fit, data.frame(genes = 20000)), "s\n")
