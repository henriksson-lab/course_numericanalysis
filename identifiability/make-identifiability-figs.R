# Figures for the identifiability/ unit ("Identifiability and uncertainty of fitted parameters").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd identifiability && Rscript make-identifiability-figs.R
#
# Base R graphics only (grDevices::svg) plus numDeriv for Jacobians.
# Colours are Okabe-Ito (colour-blind safe). Makes:
#   fig01-same-fit.svg      sparse design (first sample at 3 h): three ka values fit equally well
#   fig02-flipflop.svg      the ka/ke flip-flop: identical curves, and a ka profile with two wells
#   fig03-sse-slice.svg     SSE over (ka, V), ke re-fitted: rich vs sparse design, with the Wald ellipse
#   fig04-profiles.svg      profile of ka (others re-optimised) for both designs, Wald parabola
#   fig05-sensitivity.svg   relative sensitivities p dC/dp against time, with sample times
#   fig06-sloppy.svg        eigenvalues of J^T J (log parameters): sloppy spectra

library(numDeriv)

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

fig <- function(file, w = 9, h = 5.2, expr) {
  svg(file, width = w, height = h)
  par(mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.4, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3)
  expr
  dev.off()
}

## ---- the running example (same code as on the slides) ----------------------

D <- 100                                                         # oral dose, mg
bateman <- function(t, ka, ke, V) D * ka / (V * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
# data from *Nonlinear least squares and ODE fitting*: ka = 1 /h, ke = 0.2 /h,
# V = 20 L, D = 100 mg, noise N(0, 0.1^2) mg/L, rounded to 2 decimals
tt <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12)                         # hours
y <- c(1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61)  # mg/L
late <- tt >= 3                                                  # sparse design
tl <- tt[late]; yl <- y[late]                                     # late samples

fit  <- nls(y  ~ bateman(tt, ka, ke, V), start = list(ka = 0.5, ke = 0.1, V = 10))
fits <- nls(yl ~ bateman(tl, ka, ke, V), start = list(ka = 1, ke = 0.2, V = 20))
s2  <- deviance(fit)  / (length(y)  - 3)
s2s <- deviance(fits) / (length(yl) - 3)

# profile: fix ka, re-optimise ke and V; start from a given guess
prof1 <- function(ka, ts, ys, start = list(ke = 0.2, V = 20)) {
  f <- tryCatch(nls(ys ~ bateman(ts, ka, ke, V), start = start), error = function(e) NULL)
  if (is.null(f)) c(NA, NA, NA) else c(deviance(f), coef(f))
}
tfine <- seq(0, 12, by = 0.02)

## ---- fig01: three ka values, same fit to the late data ----------------------

fig("fig01-same-fit.svg", w = 10, h = 5.2, expr = {
  kas <- c(coef(fits)["ka"], 3, 20)
  cols <- c(oi["blue"], oi["orange"], oi["purple"])
  plot(NA, xlim = c(0, 12), ylim = c(0, 5), xlab = "time after dose (h)",
       ylab = "concentration (mg/L)")
  for (i in 1:3) {
    pr <- prof1(kas[i], tl, yl)
    lines(tfine, bateman(tfine, kas[i], pr[2], pr[3]), col = cols[i], lwd = 3, lty = i)
  }
  points(tt[!late], y[!late], pch = 1, cex = 1.6, col = oi["grey"], lwd = 2)
  points(tl, yl, pch = 19, cex = 1.6)
  abline(v = 3, col = oi["grey"], lty = 3)
  legend("topright", bty = "n", cex = 1.25, lwd = 3, lty = 1:3, col = cols,
         legend = c(sprintf("ka = %.2f /h (best fit)", kas[1]),
                    "ka = 3 /h", "ka = 20 /h"))
  legend("right", bty = "n", cex = 1.25, pch = c(19, 1), col = c("black", oi["grey"]),
         pt.cex = 1.5, legend = c("measured (from 3 h)", "not measured"))
})

## ---- fig02: flip-flop --------------------------------------------------------

fig("fig02-flipflop.svg", w = 12, h = 5.2, expr = {
  par(mfrow = c(1, 2))
  p1 <- coef(fit)
  p2 <- c(p1[2], p1[1], p1[3] * p1[2] / p1[1])
  plot(tt, y, pch = 19, cex = 1.5, xlim = c(0, 12), ylim = c(0, 4),
       xlab = "time after dose (h)", ylab = "concentration (mg/L)")
  lines(tfine, bateman(tfine, p1[1], p1[2], p1[3]), col = oi["blue"], lwd = 5)
  lines(tfine, bateman(tfine, p2[1], p2[2], p2[3]), col = oi["orange"], lwd = 2.5, lty = 2)
  legend("topright", bty = "n", cex = 1.15, lwd = c(5, 2.5), lty = 1:2,
         col = c(oi["blue"], oi["orange"]),
         legend = c(sprintf("ka %.3f, ke %.3f, V %.1f", p1[1], p1[2], p1[3]),
                    sprintf("ka %.3f, ke %.3f, V %.2f", p2[1], p2[2], p2[3])))
  # profile of ka over both branches
  g <- exp(seq(log(0.08), log(30), length.out = 121))
  g <- g[abs(g - p1[2]) > 0.004]
  d <- sapply(g, function(k) {
    a <- prof1(k, tt, y, list(ke = 0.2, V = 20))[1]
    b <- prof1(k, tt, y, list(ke = 1, V = 4))[1]
    min(a, b, na.rm = TRUE)
  })
  plot(g, (d - deviance(fit)) / s2, type = "l", lwd = 3, col = oi["blue"], log = "x",
       ylim = c(0, 100), xlab = "ka (1/h), fixed", ylab = expression(Delta * SSE / s^2))
  abline(h = 3.84, lty = 2, col = oi["grey"], lwd = 2)
  abline(v = c(p1[1], p1[2]), lty = 3, col = oi["grey"])
  text(p1[2], 92, expression(ka == hat(k)[e]), pos = 4, cex = 1.2)
  text(p1[1], 92, expression(ka == hat(k)[a]), pos = 4, cex = 1.2)
  text(25, 8, "3.84", cex = 1.15, col = "grey40")
})

## ---- fig03: SSE over (ka, V), ke re-fitted ----------------------------------

slice <- function(ts, ys, f, s2, title) {
  p <- coef(f)
  lka <- seq(log(0.2), log(30), length.out = 120)
  lV  <- seq(log(8), log(40), length.out = 120)
  # at each (ka, V) minimise the SSE over ke (kept below ka: no flip-flop branch)
  z <- outer(lka, lV, Vectorize(function(a, v)
    optimize(function(lke) sum((ys - bateman(ts, exp(a), exp(lke), exp(v)))^2),
             c(log(0.01), a - 1e-3))$objective))
  z <- (z - deviance(f)) / s2
  contour(lka, lV, z, levels = c(3.84, 10, 30, 100), axes = FALSE,
          labcex = 1, lwd = 2, col = oi["blue"], xlab = "ka (1/h), log scale",
          ylab = "V (L), log scale", main = title, cex.main = 1.4, font.main = 1)
  axis(1, at = log(c(0.2, 0.5, 1, 2, 5, 10, 20)), labels = c(0.2, 0.5, 1, 2, 5, 10, 20))
  axis(2, at = log(c(10, 15, 20, 30, 40)), labels = c(10, 15, 20, 30, 40)); box()
  # Wald ellipse for (ka, V) from the full covariance matrix vcov(f)
  C <- vcov(f)[c("ka", "V"), c("ka", "V")]
  th <- seq(0, 2 * pi, length.out = 400)
  E <- t(p[c("ka", "V")] + t(chol(C)) %*% rbind(cos(th), sin(th)) * sqrt(3.84))
  ok <- E[, 1] > 0 & E[, 2] > 0
  E[!ok, ] <- NA
  lines(log(E[, 1]), log(E[, 2]), col = oi["red"], lwd = 2.5, lty = 2)
  points(log(p["ka"]), log(p["V"]), pch = 19, cex = 1.5)
}
fig("fig03-sse-slice.svg", w = 12, h = 5.4, expr = {
  par(mfrow = c(1, 2))
  slice(tt, y, fit, s2, "all 12 samples")
  legend("bottomright", bty = "n", cex = 1.15, lwd = 2.5, lty = c(1, 2),
         col = c(oi["blue"], oi["red"]), legend = c(expression(Delta * SSE / s^2), "Wald ellipse"),
         bg = "white")
  slice(tl, yl, fits, s2s, "samples from 3 h only")
})

## ---- fig04: profiles of ka ----------------------------------------------------

fig("fig04-profiles.svg", w = 10, h = 5.4, expr = {
  g <- exp(seq(log(0.3), log(30), length.out = 201))
  dr <- sapply(g, function(k) prof1(k, tt, y)[1])
  ds <- sapply(g, function(k) prof1(k, tl, yl)[1])
  plot(g, (dr - deviance(fit)) / s2, type = "l", lwd = 3, col = oi["blue"], log = "x",
       ylim = c(0, 12), xlab = "ka (1/h), fixed; ke and V re-fitted",
       ylab = expression(Delta * SSE / s^2))
  lines(g, (ds - deviance(fits)) / s2s, lwd = 3, col = oi["orange"])
  se <- sqrt(diag(vcov(fit)))["ka"]; ses <- sqrt(diag(vcov(fits)))["ka"]
  lines(g, ((g - coef(fit)["ka"]) / se)^2, lty = 3, lwd = 2, col = oi["blue"])
  lines(g, ((g - coef(fits)["ka"]) / ses)^2, lty = 3, lwd = 2, col = oi["orange"])
  abline(h = 3.84, lty = 2, col = oi["grey"], lwd = 2)
  text(25, 4.4, "3.84", cex = 1.15, col = "grey40")
  legend("topright", bty = "n", cex = 1.2, lwd = c(3, 3, 2), lty = c(1, 1, 3),
         col = c(oi["blue"], oi["orange"], "black"),
         legend = c("profile, all 12 samples", "profile, samples from 3 h",
                    "Wald (from the SE), each design"))
})

## ---- fig05: sensitivities ----------------------------------------------------

fig("fig05-sensitivity.svg", w = 10, h = 5.2, expr = {
  p <- coef(fit)
  tf <- seq(0.01, 12, by = 0.02)
  Js <- jacobian(function(lp) bateman(tf, exp(lp[1]), exp(lp[2]), exp(lp[3])), log(p))
  par(mar = c(4.8, 5.6, 1.5, 1))
  matplot(tf, Js, type = "l", lty = c(1, 2, 4), lwd = 3,
          col = c(oi["orange"], oi["blue"], oi["green"]),
          xlab = "time after dose (h)", ylab = expression(p %.% partialdiff * C / partialdiff * p ~ "(mg/L)"))
  abline(h = 0, col = oi["grey"])
  rug(tt, ticksize = 0.05, lwd = 3)
  rug(tl, ticksize = 0.1, lwd = 3, col = oi["red"])
  legend("bottomright", bty = "n", cex = 1.25, lwd = 3, lty = c(1, 2, 4),
         col = c(oi["orange"], oi["blue"], oi["green"]),
         legend = c("ka", "ke", "V"))
})

## ---- fig06: sloppy eigenvalue spectra ----------------------------------------

pk2 <- function(t, ka, k10, k12, k21, V, D = 100) {
  s <- k10 + k12 + k21; q <- sqrt(s^2 - 4 * k21 * k10)
  al <- (s + q) / 2; be <- (s - q) / 2
  ka * D / V * ((k21 - al) / ((ka - al) * (be - al)) * exp(-al * t) +
                (k21 - be) / ((ka - be) * (al - be)) * exp(-be * t) +
                (k21 - ka) / ((al - ka) * (be - ka)) * exp(-ka * t))
}
fig("fig06-sloppy.svg", w = 10, h = 5.4, expr = {
  p <- coef(fit)
  e1 <- eigen(crossprod(jacobian(function(lp) bateman(tt, exp(lp[1]), exp(lp[2]), exp(lp[3])), log(p))))$values
  e2 <- eigen(crossprod(jacobian(function(lp) bateman(tl, exp(lp[1]), exp(lp[2]), exp(lp[3])), log(p))))$values
  p2 <- c(ka = 1, k10 = 0.2, k12 = 0.3, k21 = 0.15, V = 20)
  e3 <- eigen(crossprod(jacobian(function(lp) do.call(pk2, c(list(tt), as.list(exp(lp)))), log(p2))))$values
  plot(NA, xlim = c(0.5, 3.5), ylim = c(-5, 0.3), xaxt = "n",
       xlab = "", ylab = expression(log[10] ~ (lambda / lambda[max])))
  axis(1, at = 1:3, labels = c("1 compartment\n12 samples", "1 compartment\nfrom 3 h",
                               "2 compartments\n12 samples"), padj = 0.6, cex.axis = 1.2)
  for (i in 1:3) {
    e <- list(e1, e2, e3)[[i]]
    segments(i - 0.3, log10(e / e[1]), i + 0.3, log10(e / e[1]), lwd = 4,
             col = c(oi["blue"], oi["orange"], oi["purple"])[i])
  }
})

cat("done\n")
