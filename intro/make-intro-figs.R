# Figures for the intro/ unit ("Why numerical analysis?").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd intro && Rscript make-intro-figs.R
#
# Base R graphics only (grDevices::svg). All data are simulated or computed
# here; the growth-curve data are the same 13 OD600 values as on the slides and
# in the exercise. Colours are Okabe-Ito (colour-blind safe).

col_blue   <- "#0072B2"
col_verm   <- "#D55E00"
col_green  <- "#009E73"
col_orange <- "#E69F00"
col_purple <- "#CC79A7"
col_grey   <- "grey45"

open_svg <- function(file, width = 8, height = 5) {
  svg(file, width = width, height = height)
  par(mar = c(4.6, 5.2, 1.5, 1), cex.lab = 1.5, cex.axis = 1.3, las = 1)
}

## ---- the growth-curve data (same as slides / exercise) ---------------------
# Simulated with r = 0.5 per hour, K = 1.2, N0 = 0.02, noise sd 0.02 (set.seed(3)),
# rounded to 3 decimals; the plate reader "exports minutes".

t_min <- seq(0, 1440, by = 120)
od <- c(0.001, 0.047, 0.139, 0.282, 0.581, 0.859, 1.049,
        1.161, 1.152, 1.217, 1.182, 1.176, 1.185)
logistic <- function(p, t) p[2] / (1 + (p[2] / p[3] - 1) * exp(-p[1] * t))
sse <- function(p) sum((od - logistic(p, t_min))^2)
fit_bad  <- optim(c(r = 0.5,  K = 1, N0 = 0.05), sse)
fit_good <- optim(c(r = 0.01, K = 1, N0 = 0.05), sse)
cat("bad fit: ", signif(fit_bad$par, 4),  " SSE", signif(fit_bad$value, 4),
    " convergence", fit_bad$convergence, "\n")
cat("good fit:", signif(fit_good$par, 4), " SSE", signif(fit_good$value, 4),
    " convergence", fit_good$convergence, "\n")

## ---- fig01: two "converged" fits of the same growth curve ------------------

open_svg("fig01-growth-two-fits.svg", 8, 4.8)
tt <- seq(0, 1440, length.out = 600)
plot(t_min / 60, od, pch = 19, cex = 1.6, xlab = "time (hours)",
     ylab = "OD600", ylim = c(0, 1.35), xlim = c(0, 24))
lines(tt / 60, logistic(fit_good$par, tt), col = col_blue, lwd = 3)
lines(tt / 60, logistic(fit_bad$par, tt),  col = col_verm, lwd = 3, lty = 2)
legend("bottomright", bty = "n", cex = 1.3, lwd = 3, lty = c(1, 2),
       col = c(col_blue, col_verm),
       legend = c(sprintf("start r = 0.01: doubling %.0f min, SSE %.3f",
                          log(2) / fit_good$par["r"], fit_good$value),
                  sprintf("start r = 0.5:   doubling %.0f min, SSE %.2f",
                          log(2) / fit_bad$par["r"], fit_bad$value)))
dev.off()

## ---- fig02: Euler on drug elimination, unstable vs stable step -------------

k <- 0.8; C0 <- 10
euler <- function(h, tend) {
  n <- round(tend / h); C <- numeric(n + 1); C[1] <- C0
  for (i in 1:n) C[i + 1] <- C[i] - h * k * C[i]
  list(t = (0:n) * h, C = C)
}
e3  <- euler(3, 18)
e05 <- euler(0.5, 18)
open_svg("fig02-euler-unstable.svg", 8, 4.8)
plot(NA, xlim = c(0, 18), ylim = c(-60, 45), xlab = "time (hours)",
     ylab = "drug conc. (mg/L)")
abline(h = 0, col = "grey80")
curve(C0 * exp(-k * x), 0, 18, add = TRUE, lwd = 3, col = "black")
lines(e05$t, e05$C, type = "o", pch = 16, cex = 0.7, lwd = 2, col = col_blue)
lines(e3$t, e3$C, type = "o", pch = 15, cex = 1.4, lwd = 2.5, col = col_verm)
legend("bottomleft", bty = "n", cex = 1.3, lwd = c(3, 2, 2.5),
       pch = c(NA, 16, 15), col = c("black", col_blue, col_verm),
       legend = c("exact  10 exp(-0.8 t)", "Euler, step 0.5 h  (factor 0.6)",
                  "Euler, step 3 h  (factor -1.4)"))
dev.off()

## ---- fig03: error vs cost, and the rounding floor --------------------------

svg("fig03-error-vs-cost.svg", width = 11, height = 4.8)
par(mfrow = c(1, 2), mar = c(4.6, 5.6, 2.2, 1), cex.lab = 1.5, cex.axis = 1.3,
    las = 1, cex.main = 1.4)
nterm <- 1:17
err_series <- sapply(nterm, function(n) abs(sum(1 / factorial(0:(n - 1))) - exp(1)))
plot_err <- err_series
plot(nterm, plot_err, log = "y", type = "o", pch = 19, cex = 1.2, lwd = 2,
     col = col_blue, xlab = "number of terms (cost)", ylab = "", yaxt = "n",
     ylim = c(1e-17, 1), main = "Taylor series for e")
axis(2, at = 10^c(-16, -12, -8, -4, 0),
     labels = expression(10^-16, 10^-12, 10^-8, 10^-4, 1))
mtext("|error|", side = 2, line = 4.2, cex = 1.5, las = 0)
abline(h = .Machine$double.eps, lty = 3, col = col_grey, lwd = 2)
text(18, 1e-15, "rounding level", adj = 1, cex = 1.2, col = col_grey)

h <- 10^seq(-15, -1, by = 0.25)
err_fd <- abs((exp(h) - 1) / h - 1)
plot(h, err_fd, log = "xy", type = "o", pch = 19, cex = 0.9, lwd = 2,
     col = col_verm, xlab = "step h", ylab = "", xaxt = "n", yaxt = "n",
     main = "slope (exp(h) - 1) / h  at 0")
axis(1, at = 10^c(-15, -10, -5, -1), labels = expression(10^-15, 10^-10, 10^-5, 10^-1))
axis(2, at = 10^c(-8, -6, -4, -2, 0), labels = expression(10^-8, 10^-6, 10^-4, 10^-2, 1))
mtext("|error|", side = 2, line = 4.2, cex = 1.5, las = 0)
text(1e-4, 1e-8, "truncation\n(h too big)", cex = 1.2, col = col_grey)
text(1e-13, 1e-7, "rounding\n(h too small)", cex = 1.2, col = col_grey)
dev.off()

## ---- fig04: Taylor polynomials of exp around 0 -----------------------------

open_svg("fig04-taylor-exp.svg", 8, 5)
x <- seq(-2, 2, length.out = 400)
plot(x, exp(x), type = "l", lwd = 4, xlab = "x", ylab = "", ylim = c(-1, 7.5))
abline(h = 0, v = 0, col = "grey85")
lines(x, 1 + x, lwd = 3, col = col_verm, lty = 2)
lines(x, 1 + x + x^2 / 2, lwd = 3, col = col_blue, lty = 4)
lines(x, 1 + x + x^2 / 2 + x^3 / 6, lwd = 3, col = col_green, lty = 5)
legend("topleft", bty = "n", cex = 1.3, lwd = c(4, 3, 3, 3), lty = c(1, 2, 4, 5),
       col = c("black", col_verm, col_blue, col_green),
       legend = c(expression(e^x), expression(P[1] == 1 + x),
                  expression(P[2] == 1 + x + x^2 / 2),
                  expression(P[3] == P[2] + x^3 / 6)))
dev.off()

## ---- fig05: Taylor error vs x on log-log axes ------------------------------

open_svg("fig05-taylor-error-loglog.svg", 8, 5)
par(mar = c(4.6, 5.8, 1.5, 1))
xs <- 10^seq(-2, 0, length.out = 30)
P <- function(n, x) sapply(x, function(xx) sum(xx^(0:n) / factorial(0:n)))
cols <- c(col_verm, col_blue, col_green)
plot(NA, xlim = range(xs), ylim = c(1e-10, 1), log = "xy", xlab = "x",
     ylab = "", xaxt = "n", yaxt = "n")
axis(1, at = c(0.01, 0.03, 0.1, 0.3, 1), labels = c("0.01", "0.03", "0.1", "0.3", "1"))
axis(2, at = 10^c(-10, -8, -6, -4, -2, 0),
     labels = expression(10^-10, 10^-8, 10^-6, 10^-4, 10^-2, 1))
mtext("|error|", side = 2, line = 4.4, cex = 1.5, las = 0)
for (n in 1:3) {
  lines(xs, abs(exp(xs) - P(n, xs)), lwd = 3, col = cols[n])
  text(0.012, abs(exp(0.012) - P(n, 0.012)) * 10, adj = 0, cex = 1.3, col = cols[n],
       labels = bquote(P[.(n)] * ":  slope " * .(n + 1)))
}
dev.off()

## ---- fig06: the SSE landscape in r is flat where the bad start sits ---------

open_svg("fig06-sse-plateau.svg", 8, 4.8)
rs <- 10^seq(-3.3, 0.3, length.out = 300)
prof <- sapply(rs, function(r) sse(c(r, fit_good$par["K"], fit_good$par["N0"])))
plot(rs, prof, log = "x", type = "l", lwd = 3, col = "black",
     xlab = "growth rate r (per minute)", ylab = "SSE", ylim = c(0, 7), xaxt = "n")
axis(1, at = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1"))
points(fit_good$par["r"], fit_good$value, pch = 19, cex = 2, col = col_blue)
points(0.5, sse(c(0.5, fit_good$par["K"], fit_good$par["N0"])), pch = 17,
       cex = 2, col = col_verm)
text(0.0058, 0.25, "true minimum", adj = 1, col = col_blue, cex = 1.3)
text(0.5, 2.4, "start r = 0.5:\nflat, no slope\nto follow", col = col_verm, cex = 1.3)
dev.off()
