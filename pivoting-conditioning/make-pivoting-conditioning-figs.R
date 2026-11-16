# Figures for the pivoting-conditioning/ unit ("Pivoting, norms and conditioning").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd pivoting-conditioning && Rscript make-pivoting-conditioning-figs.R
#
# Base R only (grDevices::svg needs cairo). All data are simulated or exact
# formulas; the code mirrors the code on the slides and in the exercise.
# Writes fig01-...svg to fig07-...svg into this folder.

## ---- colours (Okabe-Ito, colour-blind safe) ----------------------------------
col_bad  <- "#D55E00"   # vermillion: no pivoting / ill-conditioned
col_good <- "#0072B2"   # blue: pivoting / well-conditioned
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
rel <- function(xhat, x) sqrt(sum((xhat - x)^2)) / sqrt(sum(x^2))

# Elimination with and without partial pivoting (as on the slides).
solve_ge <- function(A, b, pivot = TRUE) {
  n <- nrow(A)
  for (k in 1:(n - 1)) {
    if (pivot) {
      p <- which.max(abs(A[k:n, k])) + k - 1
      A[c(k, p), ] <- A[c(p, k), ]; b[c(k, p)] <- b[c(p, k)]
    }
    for (i in (k + 1):n) {
      m <- A[i, k] / A[k, k]
      A[i, ] <- A[i, ] - m * A[k, ]
      b[i]   <- b[i]   - m * b[k]
    }
  }
  x <- numeric(n)
  for (i in n:1) x[i] <- (b[i] - sum(A[i, -(1:i)] * x[-(1:i)])) / A[i, i]
  x
}

# Gaussian emission spectra: detector i sits at the peak of dye i.
panel <- function(n, spacing, width = 30) {
  peak <- 500 + spacing * (0:(n - 1))
  outer(peak, peak, function(d, p) exp(-(d - p)^2 / (2 * width^2)))
}

## ---- fig01: two dyes, nearly parallel lines ----------------------------------
open_svg("fig01-two-lines.svg", w = 7, h = 5.5)
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
       lwd = c(3, 2), lty = c(1, 2), legend = c("detector 1 (blue), 2 (orange)", "b off by 0.5 %"))
invisible(dev.off())

## ---- fig02: emission spectra and the mixing matrix ---------------------------
wl <- seq(420, 640, by = 1)
spec <- function(peak, width = 30) exp(-(wl - peak)^2 / (2 * width^2))
open_svg("fig02-spectra.svg", w = 12, h = 4.8)
par(mfrow = c(1, 2))
show_pair <- function(d2, main) {
  plot(wl, spec(500), type = "n", ylim = c(0, 1.12),
       xlab = "wavelength (nm)", ylab = "emission (relative)", main = main)
  polygon(c(wl, rev(wl)), c(spec(500), 0 * wl), col = adjustcolor(col_good, 0.25), border = NA)
  polygon(c(wl, rev(wl)), c(spec(500 + d2), 0 * wl), col = adjustcolor(col_4, 0.25), border = NA)
  lines(wl, spec(500), col = col_good, lwd = 3)
  lines(wl, spec(500 + d2), col = col_4, lwd = 3, lty = 2)
  abline(v = c(500, 500 + d2), lty = 3, lwd = 2, col = col_grey)
  ov <- exp(-d2^2 / (2 * 30^2))
  text(600, 0.95, sprintf("overlap c = %.2f", ov), cex = 1.3)
  text(600, 0.80, sprintf("kappa = %.0f", (1 + ov) / (1 - ov)), cex = 1.3)
  legend("topleft", bty = "n", cex = 1.1, lwd = c(3, 3, 2), lty = c(1, 2, 3),
         col = c(col_good, col_4, col_grey), legend = c("dye 1", "dye 2", "detectors"))
}
show_pair(30 * sqrt(-2 * log(0.5)), "peaks 35 nm apart")
show_pair(30 * sqrt(-2 * log(0.95)), "peaks 10 nm apart")
invisible(dev.off())

## ---- fig03: error without and with pivoting ----------------------------------
delta <- 10^seq(-1, -20, by = -0.25)
err_np <- err_p <- numeric(length(delta))
for (i in seq_along(delta)) {
  A <- matrix(c(delta[i], 1, 1, 1), 2, byrow = TRUE)
  b <- c(1, 2)
  x <- c(1 / (1 - delta[i]), (1 - 2 * delta[i]) / (1 - delta[i]))   # exact
  err_np[i] <- rel(solve_ge(A, b, pivot = FALSE), x)
  err_p[i]  <- rel(solve_ge(A, b, pivot = TRUE), x)
}
err_np[err_np == 0] <- 1e-17; err_p[err_p == 0] <- 1e-17
open_svg("fig03-pivot-error.svg")
plot(delta, err_np, log = "xy", type = "n", xaxt = "n", yaxt = "n",
     ylim = c(1e-17, 10), xlim = c(1e-20, 0.1),
     xlab = expression("pivot " * delta * "   (A = [" * delta * " 1; 1 1])"),
     ylab = "relative error of x")
log_axis(1, 10^seq(-20, 0, by = 4)); log_axis(2, 10^seq(-16, 0, by = 4))
abline(h = 10^seq(-16, 0, by = 4), col = "grey90")
lines(delta, pmin(eps / delta, 1), lty = 2, col = col_grey, lwd = 2)
points(delta, err_np, pch = 17, col = col_bad, cex = 1.2)
points(delta, err_p, pch = 16, col = col_good, cex = 1.2)
legend("topright", bty = "n", cex = 1.2, pch = c(17, 16, NA), lty = c(NA, NA, 2),
       lwd = c(NA, NA, 2), col = c(col_bad, col_good, col_grey),
       legend = c("no pivoting", "partial pivoting", expression(epsilon / delta)))
mtext("points at the bottom: exact result", side = 3, line = 0.5, cex = 1.2)
invisible(dev.off())

## ---- fig04: unit balls of the 1-, 2- and inf-norms ---------------------------
open_svg("fig04-unit-balls.svg", w = 6, h = 6)
par(mar = c(4.5, 5, 1, 1), mgp = c(3, 1, 0))
plot(NA, xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.75), asp = 1,
     xlab = expression(v[1]), ylab = expression(v[2]))
abline(h = 0, v = 0, col = "grey85")
th <- seq(0, 2 * pi, length.out = 300)
polygon(c(1, 0, -1, 0), c(0, 1, 0, -1), border = col_bad, lwd = 3)
lines(cos(th), sin(th), col = col_good, lwd = 3, lty = 2)
polygon(c(1, -1, -1, 1), c(1, 1, -1, -1), border = col_3, lwd = 3, lty = 4)
legend("top", bty = "n", cex = 1.2, lwd = 3, lty = c(1, 2, 4), horiz = TRUE,
       col = c(col_bad, col_good, col_3), x.intersp = 0.5,
       legend = expression(group("||", v, "||")[1] == 1, group("||", v, "||")[2] == 1,
                           group("||", v, "||")[infinity] == 1))
invisible(dev.off())

## ---- fig05: A maps the unit circle to an ellipse -----------------------------
open_svg("fig05-ellipse.svg", w = 12, h = 5.2)
par(mfrow = c(1, 2), mar = c(4.5, 5.5, 3, 1))
th <- seq(0, 2 * pi, length.out = 400)
circ <- rbind(cos(th), sin(th))
show_map <- function(c, main) {
  A <- matrix(c(1, c, c, 1), 2, byrow = TRUE)
  E <- A %*% circ
  plot(NA, xlim = c(-2.1, 2.1), ylim = c(-2.1, 2.1), asp = 1,
       xlab = expression((A * v)[1]), ylab = expression((A * v)[2]), main = main)
  abline(h = 0, v = 0, col = "grey85")
  lines(circ[1, ], circ[2, ], col = col_grey, lwd = 2, lty = 2)
  polygon(E[1, ], E[2, ], border = col_bad, lwd = 3, col = adjustcolor(col_bad, 0.15))
  s1 <- 1 + c; s2 <- 1 - c
  arrows(0, 0, s1 / sqrt(2), s1 / sqrt(2), length = 0.12, lwd = 2.5, col = col_good)
  arrows(0, 0, s2 / sqrt(2), -s2 / sqrt(2), length = 0.08, lwd = 2.5, col = col_3)
  text(s1 / sqrt(2), s1 / sqrt(2) + 0.25, bquote(sigma[max] == .(s1)), cex = 1.3, col = col_good)
  text(0.9, -1.0, bquote(sigma[min] == .(s2)), cex = 1.3, col = col_3)
  text(-1.2, 1.7, bquote(kappa == .(round(s1 / s2))), cex = 1.4)
}
show_map(0.5, "A = [1 0.5; 0.5 1]")
show_map(0.95, "A = [1 0.95; 0.95 1]")
invisible(dev.off())

## ---- fig06: digits lost vs condition number ----------------------------------
set.seed(1)
spacings <- c(60, 50, 40, 35, 30, 25, 20, 17, 15, 12, 10, 8, 7, 6, 5, 4.5, 4)
xtrue <- c(5, 1, 3, 2, 4, 1)
kap <- err_round <- err_noise <- numeric(length(spacings))
for (i in seq_along(spacings)) {
  A <- panel(6, spacings[i])
  kap[i] <- kappa(A, exact = TRUE)
  b <- drop(A %*% xtrue)
  err_round[i] <- rel(solve(A, b), xtrue)
  err_noise[i] <- median(replicate(200,
    rel(solve(A, b * (1 + rnorm(6, sd = 0.01))), xtrue)))
}
open_svg("fig06-digits-lost.svg", w = 9, h = 5.2)
par(mar = c(5, 5.5, 1, 1))
plot(kap, err_round, log = "xy", type = "n", xaxt = "n", yaxt = "n",
     ylim = c(1e-17, 1e3), xlim = c(1, 1e12),
     xlab = expression(kappa(A) ~ "  (6 dyes, peaks closer together to the right)"),
     ylab = "relative error of unmixed x")
log_axis(1, 10^seq(0, 12, by = 2)); log_axis(2, 10^seq(-16, 2, by = 3))
abline(h = 10^seq(-16, 2, by = 3), col = "grey90")
lines(c(1, 1e14), eps * c(1, 1e14), lty = 2, lwd = 2, col = col_good)
lines(c(1, 1e14), 0.01 * c(1, 1e14), lty = 2, lwd = 2, col = col_bad)
points(kap, pmax(err_round, 1e-17), pch = 16, cex = 1.3, col = col_good)
points(kap, err_noise, pch = 17, cex = 1.3, col = col_bad)
legend("bottomright", bty = "n", cex = 1.15, pch = c(17, NA, 16, NA), lty = c(NA, 2, NA, 2),
       lwd = 2, col = c(col_bad, col_bad, col_good, col_good),
       legend = c("1 % detector noise (median)", expression(kappa %*% 0.01),
                  "rounding only (exact b)", expression(kappa %*% epsilon)))
invisible(dev.off())

## ---- fig07: noise amplification in unmixed cells -----------------------------
set.seed(2)
noise <- matrix(rnorm(2 * 1000, sd = 0.01), 2)      # 1 % detector noise, 1000 cells
open_svg("fig07-unmix-cells.svg", w = 12, h = 5.2)
par(mfrow = c(1, 2), mar = c(4.8, 5.5, 3, 1))
show_cells <- function(c) {
  A <- matrix(c(1, c, c, 1), 2, byrow = TRUE)
  X <- solve(A, drop(A %*% c(1, 1)) + noise)
  plot(X[1, ], X[2, ], pch = 16, cex = 0.6, col = adjustcolor(col_good, 0.5),
       xlim = c(0.5, 1.5), ylim = c(0.5, 1.5), asp = 1,
       xlab = "unmixed dye 1", ylab = "unmixed dye 2",
       main = sprintf("overlap %.2f,  kappa = %.0f", c, (1 + c) / (1 - c)))
  points(1, 1, pch = 3, cex = 2.5, lwd = 3, col = col_bad)
  text(0.5, 1.45, sprintf("sd = %.3f", sd(X[1, ])), pos = 4, cex = 1.3)
}
show_cells(0.5)
show_cells(0.95)
invisible(dev.off())

cat("wrote fig01-fig07\n")
