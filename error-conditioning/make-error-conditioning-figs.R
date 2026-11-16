# Figures for the error-conditioning/ unit ("Error, conditioning and stability").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd error-conditioning && Rscript make-error-conditioning-figs.R
#
# Base R only (grDevices::svg needs cairo). All data are simulated or exact
# formulas; the code mirrors the code on the slides and in the exercise.
# Writes fig01-...svg to fig05-...svg into this folder.

## ---- colours (Okabe-Ito, colour-blind safe) ----------------------------------
col_bad  <- "#D55E00"   # vermillion: unstable / naive
col_good <- "#0072B2"   # blue: stable / rewritten
col_3    <- "#009E73"   # green
col_4    <- "#E69F00"   # orange
col_grey <- "grey40"

open_svg <- function(file, w = 8, h = 5) {
  svg(file, width = w, height = h)
  par(mar = c(5, 5.5, 2.5, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5,
      las = 1, mgp = c(3.6, 1, 0))
}
log_axis <- function(side, at) {
  axis(side, at = at, labels = parse(text = paste0("10^", log10(at))))
}
eps <- .Machine$double.eps

## ---- fig01: digits lost by the one-pass variance ----------------------------
# Data: mean m, sd 1 (so m/sd = m). Correct digits = -log10(relative error).
set.seed(1)
z <- rnorm(1000)
z <- (z - mean(z)) / sd(z)               # exactly sd 1 (up to rounding)
ratios <- 10^seq(0, 9, by = 0.25)
digits <- function(est, true) pmin(16, -log10(pmax(abs(est - true) / true, 1e-17)))
d_one <- d_two <- numeric(length(ratios))
for (i in seq_along(ratios)) {
  x <- ratios[i] + z
  n <- length(x)
  zs <- x - ratios[i]                     # exact: the data as actually stored
  true <- sum((zs - mean(zs))^2) / n      # shift-free reference
  d_one[i] <- digits(sum(x^2) / n - mean(x)^2, true)
  d_two[i] <- digits(sum((x - mean(x))^2) / n, true)
}
open_svg("fig01-onepass-digits.svg")
plot(ratios, d_one, log = "x", type = "n", ylim = c(0, 16.5), xaxt = "n",
     xlab = "mean / standard deviation of the data",
     ylab = "correct digits in variance")
log_axis(1, 10^seq(0, 9, by = 1))
abline(h = seq(0, 16, by = 4), col = "grey90")
lines(ratios, pmax(0, 16 - 2 * log10(ratios)), lty = 2, col = col_grey, lwd = 2)
lines(ratios, 16 - log10(2 * ratios), lty = 3, col = col_3, lwd = 3)   # conditioning limit
lines(ratios, d_two, type = "o", pch = 16, col = col_good, lwd = 2.5)
lines(ratios, d_one, type = "o", pch = 17, col = col_bad, lwd = 2.5)
legend("bottomleft", bty = "n", cex = 1.2, lwd = c(2.5, 2.5, 3, 2),
       lty = c(1, 1, 3, 2), pch = c(16, 17, NA, NA),
       col = c(col_good, col_bad, col_3, col_grey),
       legend = c("two-pass: sum((x - mean)^2)/n", "one-pass: sum(x^2)/n - mean^2",
                  "conditioning limit: 16 - log10(2 mean/sd)",
                  "16 - 2 log10(mean/sd)"))
invisible(dev.off())

## ---- fig02: cancellation in (1 + d) - 1 --------------------------------------
d <- 10^seq(-16, 0, by = 0.125)
relerr <- abs(((1 + d) - 1) - d) / d
relerr[relerr == 0] <- 1e-17
open_svg("fig02-cancellation.svg")
plot(d, relerr, log = "xy", type = "n", xaxt = "n", yaxt = "n",
     ylim = c(1e-17, 10), xlab = "d   (computing (1 + d) - 1)",
     ylab = "relative error of result")
log_axis(1, 10^seq(-16, 0, by = 2)); log_axis(2, 10^seq(-16, 0, by = 4))
abline(h = 10^seq(-16, 0, by = 4), col = "grey90")
lines(d, eps / d, lty = 2, col = col_grey, lwd = 2)
points(d, relerr, pch = 16, col = col_bad, cex = 1.1)
text(1e-7, 3e-7, expression(epsilon / d), cex = 1.5, col = col_grey, pos = 4)
mtext("points at the bottom: exact result", side = 3, line = 0.5, cex = 1.2)
invisible(dev.off())

## ---- fig03: condition numbers of simple functions ----------------------------
xx <- seq(0.02, 3, length.out = 600)
k_sub  <- abs(xx / (xx - 1))
k_log  <- abs(1 / log(xx))
k_sqrt <- rep(0.5, length(xx))
k_exp  <- abs(xx)
open_svg("fig03-kappa.svg")
plot(xx, k_sub, log = "y", type = "n", ylim = c(0.1, 1e4), yaxt = "n",
     xlab = "x", ylab = expression(kappa(x) == group("|", x * f * "'" * (x) / f(x), "|")))
log_axis(2, 10^(-1:4))
abline(h = 1, col = "grey80"); abline(v = 1, col = "grey80", lty = 3)
lines(xx, k_sub,  col = col_bad,  lwd = 3)
lines(xx, k_log,  col = col_4,    lwd = 3, lty = 2)
lines(xx, k_sqrt, col = col_good, lwd = 3)
lines(xx, k_exp,  col = col_3,    lwd = 3, lty = 4)
legend("topright", bty = "n", cex = 1.3, lwd = 3, lty = c(1, 2, 1, 4),
       col = c(col_bad, col_4, col_good, col_3),
       legend = c("f(x) = x - 1", "f(x) = log x", "f(x) = sqrt(x)", "f(x) = exp(x)"))
invisible(dev.off())

## ---- fig04: inferring dose from a response near the plateau -----------------
K <- 10                                   # half-saturation dose (uM)
resp <- function(D) D / (K + D)
dose <- function(R) K * R / (1 - R)
D <- 10^seq(-1, 4, length.out = 400)
open_svg("fig04-dose-plateau.svg", w = 9, h = 5)
plot(D, resp(D), log = "x", type = "l", lwd = 3, xaxt = "n", ylim = c(0, 1.02),
     xlab = "dose (uM)", ylab = "response (fraction of max)")
log_axis(1, 10^(-1:4))
band <- function(R0, dR, col) {
  Dlo <- dose(R0 - dR); Dhi <- dose(R0 + dR)
  y0 <- par("usr")[3]; x0 <- 10^par("usr")[1]
  rect(x0, R0 - dR, Dhi, R0 + dR, col = adjustcolor(col, 0.3), border = NA)
  rect(Dlo, y0, Dhi, R0 + dR, col = adjustcolor(col, 0.3), border = NA)
  segments(Dlo, R0 - dR, Dlo, y0, col = col, lwd = 2)
  segments(Dhi, R0 + dR, Dhi, y0, col = col, lwd = 2)
}
band(0.50, 0.005, col_good)
band(0.99, 0.005, col_bad)
lines(D, resp(D), lwd = 3)
text(0.12, 0.62, "R = 0.50 +/- 0.005:\ndose 9.8 to 10.2", pos = 4, cex = 1.2, col = col_good)
text(0.12, 0.86, "R = 0.99 +/- 0.005:\ndose 657 to 1990", pos = 4, cex = 1.2, col = col_bad)
invisible(dev.off())

## ---- fig05: 2x2 teaser, nearly parallel lines --------------------------------
# A x = b with A = [1 0.95; 0.95 1]: two nearly parallel lines.
open_svg("fig05-two-lines.svg", w = 7, h = 5.5)
par(mar = c(5, 5.5, 3.5, 1))
plot(NA, xlim = c(0.4, 1.6), ylim = c(0.4, 1.6), asp = 1,
     xlab = expression(x[1] ~ "(dye 1)"), ylab = expression(x[2] ~ "(dye 2)"))
x1 <- seq(0.4, 1.6, length.out = 50)
u <- par("usr"); clip(u[1], u[2], u[3], u[4])
line_y <- function(a, c, b) (b - a * x1) / c          # a x1 + c x2 = b
lines(x1, line_y(1, 0.95, 1.95), col = col_good, lwd = 3)
lines(x1, line_y(0.95, 1, 1.95), col = col_4, lwd = 3)
lines(x1, line_y(1, 0.95, 1.96), col = col_good, lwd = 2, lty = 2)
lines(x1, line_y(0.95, 1, 1.94), col = col_4, lwd = 2, lty = 2)
points(1, 1, pch = 16, cex = 2); points(1.2, 0.8, pch = 17, cex = 2, col = col_bad)
text(0.75, 0.7, "b = (1.95, 1.95)\nx = (1, 1)", cex = 1.2)
text(1.3, 1.1, "b = (1.96, 1.94)\nx = (1.2, 0.8)", cex = 1.2, col = col_bad)
arrows(0.8, 0.8, 0.97, 0.97, length = 0.1, lwd = 1.5)
arrows(1.3, 1.0, 1.21, 0.84, length = 0.1, lwd = 1.5, col = col_bad)
u <- par("usr"); clip(u[1], u[2], u[3] - 1, u[4] + 1)
legend("top", inset = c(0, -0.13), xpd = NA, horiz = TRUE, bty = "n", cex = 1.2,
       lwd = c(3, 2), lty = c(1, 2), legend = c("original b", "b changed by 0.5%"))
invisible(dev.off())

cat("wrote fig01-fig05\n")
