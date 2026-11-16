# Figures for the lu/ unit ("Gaussian elimination as LU").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd lu && OPENBLAS_NUM_THREADS=1 Rscript make-lu-figs.R
#
# OPENBLAS_NUM_THREADS=1 makes the timings single-threaded, so the slopes
# show the operation counts and not how many cores join in (see the
# "When timing misleads" slide in convergence-cost). Timings are machine
# dependent; the slopes are not.
#
# Base R graphics only (grDevices::svg) plus the Matrix package (lu()).
# Colours are Okabe-Ito (colour-blind safe). Makes:
#   fig01-mixing.svg      isotope mixing polygon: 3 diet sources, 4 consumers
#   fig02-lu-timing.svg   time of lu(A) and of one solve with the factors vs n
#   fig03-many-rhs.svg    total time for m right-hand sides: three strategies
#   fig04-inverse.svg     relative residual of solve() vs solve(A) %*% b

library(Matrix)

## Minimum over `reps` timings of f() (the minimum is least affected by other
## jobs on the machine), divided by `inner` repetitions inside one timing.
tmin <- function(f, reps = 5, inner = 1)
  min(replicate(reps, system.time(for (k in seq_len(inner)) f())["elapsed"])) / inner

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- fig01: the mixing polygon ---------------------------------------------
## Same numbers as on the slides: source signatures (corrected for trophic
## discrimination) and four consumers.

src <- data.frame(name = c("C3 plants", "maize (C4)", "fish"),
                  c13 = c(-26, -10, -18), n15 = c(2, 6, 14))
cons <- data.frame(c13 = c(-19.6, -14.8, -17.2, -24.0),
                   n15 = c(5.6, 6.8, 11.2, 9.0))

svg("fig01-mixing.svg", width = 8, height = 6)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.4, cex.axis = 1.3)
plot(NA, xlim = c(-28, -8), ylim = c(0, 16), las = 1,
     xlab = expression(delta^13 * C ~ "(‰)"),
     ylab = expression(delta^15 * N ~ "(‰)"))
polygon(src$c13, src$n15, col = adjustcolor(oi["sky"], 0.25), border = oi["blue"],
        lwd = 2)
points(src$c13, src$n15, pch = 15, cex = 2.2, col = oi["blue"])
text(src$c13, src$n15, src$name, pos = c(1, 1, 3), cex = 1.3)
inside <- c(TRUE, TRUE, TRUE, FALSE)
points(cons$c13, cons$n15, pch = 19, cex = 1.8,
       col = ifelse(inside, oi["orange"], oi["red"]))
text(cons$c13, cons$n15, paste0("#", 1:4), pos = 4, cex = 1.3)
text(-24, 9, "outside: a fraction < 0", pos = 3, cex = 1.2, col = oi["red"])
legend("topright", c("diet source", "consumer"), pch = c(15, 19),
       col = c(oi["blue"], oi["orange"]), pt.cex = c(2, 1.6), bty = "n", cex = 1.3)
dev.off()

## ---- fig02: cost of factorising vs cost of substitution --------------------

set.seed(1)
ns <- c(250, 500, 1000, 2000, 4000, 8000)
t_fac <- t_sol <- numeric(length(ns))
for (i in seq_along(ns)) {
  n <- ns[i]
  A <- matrix(rnorm(n * n), n); b <- rnorm(n)
  F <- lu(A)
  t_fac[i] <- tmin(function() lu(A), reps = if (n > 4000) 2 else 5,
                  inner = max(1, round(2000 / n)))
  t_sol[i] <- tmin(function() solve(F, b), inner = 100)
}
print(data.frame(n = ns, t_fac, t_sol))
s_fac <- coef(lm(log10(t_fac[ns >= 1000]) ~ log10(ns[ns >= 1000])))[2]
s_sol <- coef(lm(log10(t_sol[ns >= 1000]) ~ log10(ns[ns >= 1000])))[2]
cat("slopes:", s_fac, s_sol, "\n")

svg("fig02-lu-timing.svg", width = 9, height = 5.5)
par(mar = c(5, 6.5, 1, 1), cex.lab = 1.4, cex.axis = 1.3)
plot(ns, t_fac, log = "xy", type = "b", pch = 19, lwd = 3, cex = 1.5,
     col = oi["red"], ylim = range(c(t_fac, t_sol)), axes = FALSE,
     xlab = "matrix size n", ylab = "")
lines(ns, t_sol, type = "b", pch = 17, lwd = 3, cex = 1.5, col = oi["blue"])
axis(1, at = ns)
axis(2, at = 10^(-5:1), labels = c("10 µs", "100 µs", "1 ms", "10 ms",
                                   "0.1 s", "1 s", "10 s"), las = 1)
title(ylab = "seconds (one core)", line = 5)
box()
legend("topleft",
       c(sprintf("factorise: lu(A)   slope %.1f", s_fac),
         sprintf("forward + back substitution   slope %.1f", s_sol)),
       pch = c(19, 17), lwd = 3, col = c(oi["red"], oi["blue"]), bty = "n",
       cex = 1.3)
dev.off()

## ---- fig03: many right-hand sides --------------------------------------------

set.seed(2)
n <- 1000
A <- matrix(rnorm(n * n), n); b <- rnorm(n)
F <- lu(A); Ai <- solve(A)
t_solve1 <- tmin(function() solve(A, b), inner = 3)
t_lu     <- tmin(function() lu(A), inner = 3)
t_sub    <- tmin(function() solve(F, b), inner = 100)
t_inv    <- tmin(function() solve(A), inner = 3)
t_mv     <- tmin(function() Ai %*% b, inner = 100)
cat("n = 1000:", t_solve1, t_lu, t_sub, t_inv, t_mv, "\n")

m <- 10^seq(0, 4, length.out = 50)
svg("fig03-many-rhs.svg", width = 9, height = 5.5)
par(mar = c(5, 6.5, 1, 1), cex.lab = 1.4, cex.axis = 1.3)
plot(m, m * t_solve1, log = "xy", type = "l", lwd = 3, col = oi["red"],
     ylim = c(0.01, 1000), axes = FALSE,
     xlab = "number of right-hand sides m  (n = 1000)", ylab = "")
lines(m, t_lu + m * t_sub, lwd = 3, col = oi["blue"])
lines(m, t_inv + m * t_mv, lwd = 3, lty = 2, col = oi["orange"])
axis(1, at = 10^(0:4), labels = c("1", "10", "100", "1000", "10 000"))
axis(2, at = 10^(-2:3), labels = c("0.01", "0.1", "1", "10", "100", "1000"), las = 1)
title(ylab = "total seconds (one core)", line = 4.5)
box()
legend("topleft",
       c("solve(A, b) for each b: factorise m times",
         "lu(A) once, then solve(F, b) for each b",
         "solve(A) once, then Ainv %*% b for each b"),
       lwd = 3, lty = c(1, 1, 2), col = c(oi["red"], oi["blue"], oi["orange"]),
       bty = "n", cex = 1.2)
dev.off()

## ---- fig04: accuracy of the inverse ------------------------------------------
## Random 50 x 50 matrices with prescribed condition number kappa
## (A = Q1 diag(s) Q2^T with singular values from 1 down to 1/kappa).

set.seed(3)
n <- 50
kap <- 10^seq(1, 15, length.out = 60)
res_lu <- res_inv <- numeric(length(kap))
for (i in seq_along(kap)) {
  Q1 <- qr.Q(qr(matrix(rnorm(n * n), n)))
  Q2 <- qr.Q(qr(matrix(rnorm(n * n), n)))
  A <- Q1 %*% diag(10^seq(0, -log10(kap[i]), length.out = n)) %*% t(Q2)
  b <- A %*% rnorm(n)
  x1 <- solve(A, b, tol = 0)
  x2 <- solve(A, tol = 0) %*% b
  rr <- function(x) max(abs(b - A %*% x)) / (max(abs(A)) * max(abs(x)))
  res_lu[i] <- rr(x1); res_inv[i] <- rr(x2)
}

svg("fig04-inverse.svg", width = 9, height = 5.5)
par(mar = c(5, 6.5, 1, 1), cex.lab = 1.4, cex.axis = 1.3)
plot(kap, res_inv, log = "xy", pch = 17, cex = 1.3, col = oi["orange"],
     ylim = c(1e-19, 1e-2), axes = FALSE,
     xlab = expression("condition number " * kappa(A)), ylab = "")
points(kap, res_lu, pch = 19, cex = 1.3, col = oi["blue"])
axis(1, at = 10^seq(1, 15, 2), labels = parse(text = paste0("10^", seq(1, 15, 2))))
axis(2, at = 10^seq(-18, -2, 4), labels = parse(text = paste0("10^", seq(-18, -2, 4))),
     las = 1)
title(ylab = "relative residual", line = 4.5)
box()
abline(h = .Machine$double.eps, lty = 3, col = oi["grey"], lwd = 2)
text(10, .Machine$double.eps, expression(epsilon), pos = 3, cex = 1.4)
legend("topleft", c("solve(A) %*% b", "solve(A, b)"), pch = c(17, 19),
       col = c(oi["orange"], oi["blue"]), bty = "n", cex = 1.3)
dev.off()
