# Figures for the quadrature/ unit ("Numerical integration").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd quadrature && Rscript make-quadrature-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Everything is computed here with the same code and
# numbers as on the slides. Makes:
#   fig01-pk-auc.svg        sampled drug concentrations, trapezoids, model
#   fig02-rules.svg         left Riemann, midpoint, trapezoid, Simpson
#   fig03-orders.svg        error vs number of intervals (log-log), AUC 0-6 h
#   fig04-smoothness.svg    order drops for sqrt(x) and 1/sqrt(x(1-x))
#   fig05-adaptive.svg      where adaptive Simpson puts its points
#   fig06-kronrod.svg       the 21 points integrate() uses on [0, 20000]
#   fig07-substitution.svg  infinite interval and singularity, before/after
#   fig08-peak-sampling.svg chromatography peak sampled every 2 sigma / 1 sigma
#   fig09-lul.svg           linear vs log-down interpolation of PK data

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")
tr <- function(col, a = 0.35) adjustcolor(col, alpha.f = a)

## ---- shared: the pharmacokinetic example ------------------------------------
## One-compartment model, oral dose 500 mg, V = 50 L, ka = 1.5 /h, ke = 0.2 /h.
model <- function(t) 500 / 50 * 1.5 / 1.3 * (exp(-0.2 * t) - exp(-1.5 * t))
t_h  <- c(0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24)
conc <- c(0, 4.99, 6.87, 7.16, 6.20, 5.16, 3.47, 2.33, 1.05, 0.09)  # = round(model(t_h), 2)
exact6 <- 32.6243601574044                                          # AUC 0-6 h

trap <- function(f, a, b, n) {
  x <- seq(a, b, length.out = n + 1);  h <- (b - a) / n
  y <- f(x);  h * (sum(y) - (y[1] + y[n + 1]) / 2)
}
left <- function(f, a, b, n) { h <- (b - a) / n; h * sum(f(a + (0:(n - 1)) * h)) }
mid  <- function(f, a, b, n) { h <- (b - a) / n; h * sum(f(a + (1:n - 0.5) * h)) }
simpson <- function(f, a, b, n) {
  x <- seq(a, b, length.out = n + 1);  h <- (b - a) / n
  w <- c(1, rep(c(4, 2), n / 2 - 1), 4, 1);  h / 3 * sum(w * f(x))
}
setpar <- function(mar = c(5, 5.5, 1, 1))
  par(mar = mar, cex.lab = 1.5, cex.axis = 1.35, mgp = c(3.2, 1, 0))

## ---- fig01: the AUC problem --------------------------------------------------
svg("fig01-pk-auc.svg", width = 10, height = 5.2)
setpar()
tt <- seq(0, 24, length.out = 500)
plot(tt, model(tt), type = "n", xlab = "time after dose (h)",
     ylab = "concentration (mg/L)", ylim = c(0, 8), las = 1, xaxt = "n")
axis(1, at = c(0, 2, 4, 6, 8, 12, 24))
for (i in 1:9)
  polygon(t_h[c(i, i, i + 1, i + 1)], c(0, conc[i], conc[i + 1], 0),
          col = tr(oi["sky"]), border = oi["blue"])
lines(tt, model(tt), lwd = 3, col = oi["red"])
points(t_h, conc, pch = 19, cex = 1.6)
legend("topright", c("blood samples", "straight lines between samples (trapezoids)",
                     "true concentration (model)"),
       pch = c(19, NA, NA), lty = c(NA, 1, 1), lwd = c(NA, 8, 3),
       col = c("black", tr(oi["sky"], 0.6), oi["red"]), bty = "n", cex = 1.3)
text(17, 2.4, "trapezoids: 51.6 mg h/L\nmodel: 49.5 mg h/L", cex = 1.35, adj = 0)
dev.off()

## ---- fig02: four rules -------------------------------------------------------
svg("fig02-rules.svg", width = 12, height = 3.8)
par(mfrow = c(1, 4), mar = c(3, 1, 3, 0.5), cex.main = 1.7, cex.axis = 1.3)
x <- seq(0, 6, length.out = 400);  n <- 4;  h <- 6 / n;  xs <- seq(0, 6, by = h)
base <- function(main) {
  plot(x, model(x), type = "n", ylim = c(0, 8), axes = FALSE, main = main,
       xlab = "", ylab = "")
  axis(1, at = xs)
}
base("left Riemann")
for (i in 1:n) rect(xs[i], 0, xs[i + 1], model(xs[i]), col = tr(oi["sky"]), border = oi["blue"])
lines(x, model(x), lwd = 3, col = oi["red"])
base("midpoint")
for (i in 1:n) rect(xs[i], 0, xs[i + 1], model(xs[i] + h / 2), col = tr(oi["sky"]), border = oi["blue"])
points(xs[-1] - h / 2, model(xs[-1] - h / 2), pch = 19)
lines(x, model(x), lwd = 3, col = oi["red"])
base("trapezoid")
for (i in 1:n) polygon(xs[c(i, i, i + 1, i + 1)], c(0, model(xs[i]), model(xs[i + 1]), 0),
                       col = tr(oi["sky"]), border = oi["blue"])
lines(x, model(x), lwd = 3, col = oi["red"])
base("Simpson")
for (i in c(1, 3)) {
  xp <- xs[i:(i + 2)];  cf <- solve(cbind(1, xp, xp^2), model(xp))
  xx <- seq(xp[1], xp[3], length.out = 100)
  polygon(c(xx, rev(xx)), c(cf[1] + cf[2] * xx + cf[3] * xx^2, 0 * xx),
          col = tr(oi["sky"]), border = oi["blue"])
}
lines(x, model(x), lwd = 3, col = oi["red"])
points(xs, model(xs), pch = 19)
dev.off()

## ---- fig03: orders on a log-log plot -----------------------------------------
n <- 2^(1:10)
E <- sapply(list(left, mid, trap, simpson),
            function(r) sapply(n, function(k) abs(r(model, 0, 6, k) - exact6)))
svg("fig03-orders.svg", width = 9, height = 5.6)
setpar(c(5, 6.5, 1, 1))
cols <- c(oi["grey"], oi["green"], oi["blue"], oi["red"])
matplot(n, E, log = "xy", type = "b", pch = 19, lty = 1, lwd = 2.5, col = cols,
        xlab = "number of intervals n  (h = 6/n)", ylab = "", axes = FALSE)
axis(1, at = n);  axis(2, at = 10^seq(-10, 2, 2), las = 1,
                       labels = parse(text = paste0("10^", seq(-10, 2, 2))))
title(ylab = "|error| of AUC 0-6 h", line = 4.8);  box()
legend("bottomleft", c("left Riemann: slope -1", "midpoint: slope -2",
                       "trapezoid: slope -2", "Simpson: slope -4"),
       col = cols, lwd = 2.5, pch = 19, bty = "n", cex = 1.3)
dev.off()

## ---- fig04: order needs smoothness -------------------------------------------
n <- 2^(2:10)
e_pk   <- sapply(n, function(k) abs(simpson(model, 0, 6, k) - exact6))
e_sq_s <- sapply(n, function(k) abs(simpson(sqrt, 0, 1, k) - 2/3))
e_sq_t <- sapply(n, function(k) abs(trap(sqrt, 0, 1, k) - 2/3))
jeff   <- function(x) 1 / sqrt(x * (1 - x))
e_jm   <- sapply(n, function(k) abs(mid(jeff, 0, 1, k) - pi))
svg("fig04-smoothness.svg", width = 9, height = 5.6)
setpar(c(5, 6.5, 1, 1))
cols <- c(oi["red"], oi["orange"], oi["blue"], oi["purple"])
matplot(n, cbind(e_pk, e_sq_s, e_sq_t, e_jm), log = "xy", type = "b", pch = 19,
        lty = 1, lwd = 2.5, col = cols, xlab = "number of intervals n",
        ylab = "", axes = FALSE, ylim = c(1e-10, 1))
axis(1, at = n);  axis(2, at = 10^seq(-10, 0, 2), las = 1,
                       labels = parse(text = paste0("10^", seq(-10, 0, 2))))
title(ylab = "|error|", line = 4.8);  box()
legend("bottomleft", c("Simpson, smooth PK curve: slope -4",
                       expression("Simpson, " * sqrt(x) * ": slope -1.5"),
                       expression("trapezoid, " * sqrt(x) * ": slope -1.5"),
                       expression("midpoint, " * 1/sqrt(x(1-x)) * ": slope -0.5")),
       col = cols, lwd = 2.5, pch = 19, bty = "n", cex = 1.25)
dev.off()

## ---- fig05: adaptive Simpson ------------------------------------------------
adapt <- function(f, a, b, fa, fm, fb, tol) {
  m <- (a + b) / 2;  fl <- f((a + m) / 2);  fr <- f((m + b) / 2)
  S1 <- (b - a) / 6 * (fa + 4 * fm + fb)
  S2 <- (b - a) / 12 * (fa + 4 * fl + 2 * fm + 4 * fr + fb)
  if (abs(S2 - S1) < 15 * tol) return(S2 + (S2 - S1) / 15)
  adapt(f, a, m, fa, fl, fm, tol / 2) + adapt(f, m, b, fm, fr, fb, tol / 2)
}
xs <- c()
fc <- function(t) { xs <<- c(xs, t); model(t) }
invisible(adapt(fc, 0, 24, fc(0), fc(12), fc(24), 1e-6))
svg("fig05-adaptive.svg", width = 10, height = 5.2)
setpar()
tt <- seq(0, 24, length.out = 500)
plot(tt, model(tt), type = "l", lwd = 3, col = oi["red"], las = 1,
     xlab = "time after dose (h)", ylab = "concentration (mg/L)", ylim = c(-0.8, 8))
segments(xs, -0.8, xs, -0.2, col = oi["blue"], lwd = 1.2)
text(24, 7.2, paste0(length(xs), " evaluations, error 9e-10\n",
                     "(uniform Simpson: 513 evaluations, error 1e-6)"),
     adj = 1, cex = 1.3)
dev.off()

## ---- fig06: the 21 Gauss-Kronrod points on [0, 20000] ------------------------
seen <- NULL
spy <- function(x) { seen <<- c(seen, x); dnorm(x) }
invisible(integrate(spy, 0, 20000))
svg("fig06-kronrod.svg", width = 12, height = 4.6)
par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1), cex.lab = 1.5, cex.axis = 1.3,
    cex.main = 1.5)
x <- seq(0, 20000, length.out = 4001)
plot(x, dnorm(x), type = "l", lwd = 3, col = oi["red"], ylim = c(0, 0.42),
     xlab = "x", ylab = "dnorm(x)", main = "what integrate() sees on [0, 20000]", las = 1)
points(seen, dnorm(seen), pch = 19, col = oi["blue"], cex = 1.4)
x <- seq(0, 60, length.out = 600)
plot(x, dnorm(x), type = "l", lwd = 3, col = oi["red"], ylim = c(0, 0.42),
     xlab = "x", ylab = "dnorm(x)", main = "zoom: [0, 60]", las = 1)
points(seen[seen < 60], dnorm(seen[seen < 60]), pch = 19, col = oi["blue"], cex = 1.6)
text(43.4, 0.06, "first point: x = 43.4\ndnorm = 0", cex = 1.3)
dev.off()

## ---- fig07: substitution -----------------------------------------------------
svg("fig07-substitution.svg", width = 12, height = 4.6)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), cex.lab = 1.5, cex.axis = 1.3,
    cex.main = 1.5)
x <- seq(0.002, 0.998, length.out = 500)
plot(x, jeff(x), type = "l", lwd = 3, col = oi["red"], ylim = c(0, 8), las = 1,
     xlab = "x  or  u", ylab = "integrand", main = expression(1/sqrt(x(1-x)) ~ "on (0, 1)"))
th <- seq(0, pi / 2, length.out = 100)
lines(th / (pi / 2), rep(2, 100) , lwd = 3, col = oi["blue"], lty = 2)
legend("top", c("before: infinite at both ends",
                expression("after x = sin"^2 * theta * ": constant 2 on (0, " * pi/2 * ")")),
       col = c(oi["red"], oi["blue"]), lwd = 3, lty = c(1, 2), bty = "n", cex = 1.15)
x <- seq(0, 10, length.out = 500);  t <- seq(0, 0.999, length.out = 500)
plot(x / 10, exp(-x), type = "l", lwd = 3, col = oi["red"], ylim = c(0, 1.7),
     las = 1, xlab = "x / 10  or  t", ylab = "integrand",
     main = expression(e^-x ~ "on [0, " * infinity * ")"))
lines(t, exp(-t / (1 - t)) / (1 - t)^2, lwd = 3, col = oi["blue"], lty = 2)
legend("topright", c(expression("before: " * e^-x * ", x up to " * infinity),
                     expression("after x = t/(1-t): t in [0, 1)")),
       col = c(oi["red"], oi["blue"]), lwd = 3, lty = c(1, 2), bty = "n", cex = 1.15)
dev.off()

## ---- fig08: chromatography peak sampling ------------------------------------
pk <- function(t) 1000 * exp(-(t - 300)^2 / (2 * 3^2))   # sigma = 3 s
A  <- 1000 * 3 * sqrt(2 * pi)
svg("fig08-peak-sampling.svg", width = 12, height = 4.4)
par(mfrow = c(1, 3), mar = c(5, 5, 3.5, 1), cex.lab = 1.5, cex.axis = 1.3,
    cex.main = 1.45)
tt <- seq(285, 315, length.out = 400)
for (cfg in list(c(6, 0), c(6, 3), c(3, 1.5))) {
  h <- cfg[1];  xs <- seq(240 + cfg[2], 360, by = h);  ys <- pk(xs)
  a <- h * (sum(ys) - (ys[1] + ys[length(ys)]) / 2)
  plot(tt, pk(tt), type = "l", lwd = 3, col = oi["red"], ylim = c(0, 1050),
       xlab = "retention time (s)", ylab = "signal", las = 1,
       main = sprintf("every %g s: area %s, height %+.0f%%", h,
                      if (abs(a / A - 1) < 1e-4) sprintf("%+.0e", a / A - 1)
                      else sprintf("%+.1f%%", 100 * (a / A - 1)),
                      100 * (max(ys) / 1000 - 1)))
  k <- xs >= 280 & xs <= 320
  polygon(c(xs[k][1], xs[k], xs[k][sum(k)]), c(0, ys[k], 0),
          col = tr(oi["sky"]), border = oi["blue"])
  points(xs[k], ys[k], pch = 19, cex = 1.5)
}
dev.off()

## ---- fig09: linear-up / log-down ---------------------------------------------
svg("fig09-lul.svg", width = 12, height = 4.6)
par(mfrow = c(1, 2), mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
tt <- seq(0, 24, length.out = 500)
loginterp <- function(t) {               # log-linear between samples after Cmax
  i <- findInterval(t, t_h, rightmost.closed = TRUE)
  y1 <- conc[i];  y2 <- conc[i + 1];  f <- (t - t_h[i]) / (t_h[i + 1] - t_h[i])
  ifelse(y2 < y1, y1 * (y2 / y1)^f, y1 + f * (y2 - y1))
}
for (lg in c("", "y")) {
  rng <- if (lg == "y") c(0.05, 10) else c(0, 8)
  t2 <- if (lg == "y") tt[tt > 0.2] else tt
  plot(t2, model(t2), type = "l", lwd = 3, col = oi["red"], log = lg, ylim = rng,
       las = 1, xlab = "time after dose (h)", ylab = "concentration (mg/L)")
  lines(t2, approx(t_h, conc, t2)$y, lwd = 2.5, col = oi["blue"], lty = 2)
  lines(t2, loginterp(t2), lwd = 2.5, col = oi["green"], lty = 3)
  points(t_h, conc, pch = 19, cex = 1.4)
  if (lg == "")
    legend("topright", c("true curve", "linear between samples (51.6)",
                         "linear up / log down (48.9)"),
           col = c(oi["red"], oi["blue"], oi["green"]), lwd = 3, lty = 1:3,
           bty = "n", cex = 1.3)
}
dev.off()
