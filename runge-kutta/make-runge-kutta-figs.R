# Figures for the runge-kutta/ unit ("Runge–Kutta and adaptive step size").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd runge-kutta && Rscript make-runge-kutta-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Makes:
#   fig01-lv-euler-rk4.svg     predator-prey orbit: Euler spirals out, RK4 closes
#   fig02-midpoint-heun.svg    one step of Euler, midpoint and Heun on the logistic curve
#   fig03-rk4-stages.svg       where RK4 evaluates its four slopes
#   fig04-order.svg            error vs h on log-log axes: slopes 1, 2, 4
#   fig05-cost.svg             error vs number of f evaluations
#   fig06-repressilator.svg    repressilator proteins and the step sizes of an adaptive solver
#   fig07-lv-invariant.svg     drift of the conserved quantity V(t) for each method
#   fig08-tolerance.svg        achieved error vs rtol; atol too large gives a floor
#   fig09-dense-output.svg     solver steps vs dense output
#   fig10-events.svg           an event located inside a step; an event missed by a big step

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- the methods (same code as on the slides) ---------------------------------

euler_step <- function(f, t, y, h) y + h * f(t, y)
midpoint_step <- function(f, t, y, h) {
  k1 <- f(t, y)
  k2 <- f(t + h / 2, y + h / 2 * k1)
  y + h * k2
}
heun_step <- function(f, t, y, h) {
  k1 <- f(t, y)
  k2 <- f(t + h, y + h * k1)
  y + h / 2 * (k1 + k2)
}
rk4_step <- function(f, t, y, h) {
  k1 <- f(t, y)
  k2 <- f(t + h / 2, y + h / 2 * k1)
  k3 <- f(t + h / 2, y + h / 2 * k2)
  k4 <- f(t + h, y + h * k3)
  y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
}
rk_solve <- function(step, f, t0, y0, h, n) {
  t <- t0 + (0:n) * h
  Y <- matrix(0, n + 1, length(y0)); Y[1, ] <- y0
  for (i in 1:n) Y[i + 1, ] <- step(f, t[i], Y[i, ], h)
  list(t = t, y = Y)
}

# Dormand-Prince 5(4), the method behind ode45 / RK45 / Dopri5, with the
# standard controller. Records accepted steps (t, y, f at the step end).
dp5 <- function(f, t0, y0, tend, rtol = 1e-6, atol = 1e-6, h = 0.01) {
  A <- list(c(), c(1/5), c(3/40, 9/40), c(44/45, -56/15, 32/9),
            c(19372/6561, -25360/2187, 64448/6561, -212/729),
            c(9017/3168, -355/33, 46732/5247, 49/176, -5103/18656),
            c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84))
  cc <- c(0, 1/5, 3/10, 4/5, 8/9, 1, 1)
  b5 <- c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0)
  b4 <- c(5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40)
  t <- t0; y <- y0; Tt <- t; Y <- y; Fv <- f(t, y); H <- c(); rej <- c()
  while (t < tend) {
    h <- min(h, tend - t)
    k <- matrix(0, 7, length(y))
    for (s in 1:7) {
      ys <- y
      if (s > 1) for (j in 1:(s - 1)) ys <- ys + h * A[[s]][j] * k[j, ]
      k[s, ] <- f(t + cc[s] * h, ys)
    }
    y5 <- y + h * colSums(b5 * k); y4 <- y + h * colSums(b4 * k)
    sc <- atol + rtol * pmax(abs(y), abs(y5))
    err <- sqrt(mean(((y5 - y4) / sc)^2))
    if (err <= 1) {
      t <- t + h; y <- y5
      Tt <- c(Tt, t); Y <- rbind(Y, y); Fv <- rbind(Fv, k[7, ]); H <- c(H, h)
    } else rej <- c(rej, t)
    h <- h * min(5, max(0.2, 0.9 * err^(-1/5)))
  }
  list(t = Tt, y = unname(Y), f = unname(Fv), h = H, rej = rej)
}

# cubic Hermite interpolation between two accepted steps (a simple dense output)
hermite <- function(tq, t0, t1, y0, y1, f0, f1) {
  h <- t1 - t0; s <- (tq - t0) / h
  (2*s^3 - 3*s^2 + 1) * y0 + (s^3 - 2*s^2 + s) * h * f0 +
    (-2*s^3 + 3*s^2) * y1 + (s^3 - s^2) * h * f1
}

## ---- running examples -------------------------------------------------------

K <- 1.2; N0 <- 0.05; r <- 0.8                 # logistic OD curve, per hour
flog <- function(t, N) r * N * (1 - N / K)
exact <- function(t) K / (1 + (K / N0 - 1) * exp(-r * t))

a <- 1; b <- 0.1; c <- 0.075; d <- 1.5          # Lotka-Volterra, per year
lv <- function(t, y) c(a * y[1] - b * y[1] * y[2], c * y[1] * y[2] - d * y[2])
V <- function(x, y) c * x - d * log(x) + b * y - a * log(y)
lv0 <- c(10, 5)

beta <- 0.3; gamma <- 0.1                       # SIR, fractions of the population
sir <- function(t, y) c(-beta * y[1] * y[2], beta * y[1] * y[2] - gamma * y[2], gamma * y[2])

alpha <- 216; alpha0 <- 0.216; bet <- 5; nH <- 2  # repressilator, dimensionless
repress <- function(t, y) {
  m <- y[1:3]; p <- y[4:6]
  prep <- p[c(3, 1, 2)]
  c(-m + alpha / (1 + prep^nH) + alpha0, -bet * (p - m))
}
rep0 <- c(1, 0, 0, 0, 0, 0)

## ---- fig01: predator-prey, Euler vs RK4 -------------------------------------

svg("fig01-lv-euler-rk4.svg", width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(5, 5.5, 2.5, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
e <- rk_solve(euler_step, lv, 0, lv0, 0.02, 1000)
k <- rk_solve(rk4_step, lv, 0, lv0, 0.1, 200)
plot(e$y, type = "l", lwd = 2, col = oi["red"], xlab = "prey x (hares)",
     ylab = "predators y (lynx)", las = 1, main = "Euler, h = 0.02 yr")
points(lv0[1], lv0[2], pch = 19, cex = 1.5)
plot(k$y, type = "l", lwd = 2, col = oi["blue"], xlab = "prey x (hares)",
     ylab = "predators y (lynx)", las = 1, xlim = range(e$y[, 1]), ylim = range(e$y[, 2]),
     main = "RK4, h = 0.1 yr")
points(lv0[1], lv0[2], pch = 19, cex = 1.5)
dev.off()

## ---- fig02: Euler, midpoint, Heun: one step --------------------------------

svg("fig02-midpoint-heun.svg", width = 11, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
h <- 2; t0 <- 1; y0 <- exact(t0)
tt <- seq(0.6, 3.4, length.out = 200)
plot(tt, exact(tt), type = "l", lwd = 4, col = oi["grey"], xlab = "time t (h)",
     ylab = "OD  N(t)", las = 1, ylim = c(0, 0.75))
k1 <- flog(t0, y0)
ye <- y0 + h * k1
ym <- y0 + h / 2 * k1; k2m <- flog(t0 + h / 2, ym)
k2h <- flog(t0 + h, ye)
segments(t0, y0, t0 + h, ye, lwd = 3, col = oi["red"])
segments(t0, y0, t0 + h / 2, ym, lwd = 2, lty = 2, col = oi["orange"])
points(t0 + h / 2, ym, pch = 1, cex = 1.8, lwd = 2, col = oi["orange"])
segments(t0, y0, t0 + h, y0 + h * k2m, lwd = 3, col = oi["orange"])
segments(t0, y0, t0 + h, y0 + h / 2 * (k1 + k2h), lwd = 3, col = oi["green"])
points(t0 + h, ye, pch = 1, cex = 1.8, lwd = 2, col = oi["green"])
points(c(t0 + h, t0 + h, t0 + h), c(ye, y0 + h * k2m, y0 + h / 2 * (k1 + k2h)),
       pch = 19, cex = 1.5, col = oi[c("red", "orange", "green")])
points(t0, y0, pch = 19, cex = 1.5)
points(t0 + h, exact(t0 + h), pch = 4, cex = 2, lwd = 3)
legend("topleft", c("exact", "Euler: slope k1 at the start",
                    "midpoint: slope at the half-step point (open circle)",
                    "Heun: average of k1 and the slope at the Euler point",
                    "x  exact value at t + h"),
       col = c(oi[c("grey", "red", "orange", "green")], "black"), lwd = c(4, 3, 3, 3, NA),
       cex = 1.2, bty = "n")
dev.off()

## ---- fig03: RK4 stages ------------------------------------------------------

svg("fig03-rk4-stages.svg", width = 11, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(tt, exact(tt), type = "l", lwd = 4, col = oi["grey"], xlab = "time t (h)",
     ylab = "OD  N(t)", las = 1, ylim = c(0, 0.75))
k1 <- flog(t0, y0); p2 <- y0 + h / 2 * k1; k2 <- flog(t0 + h / 2, p2)
p3 <- y0 + h / 2 * k2; k3 <- flog(t0 + h / 2, p3); p4 <- y0 + h * k3; k4 <- flog(t0 + h, p4)
seg <- function(tp, yp, k, col) {
  dt <- 0.25; segments(tp - dt, yp - dt * k, tp + dt, yp + dt * k, lwd = 4, col = col)
  points(tp, yp, pch = 19, cex = 1.6, col = col)
}
cols <- oi[c("red", "orange", "purple", "blue")]
seg(t0, y0, k1, cols[1]); seg(t0 + h / 2, p2, k2, cols[2])
seg(t0 + h / 2, p3, k3, cols[3]); seg(t0 + h, p4, k4, cols[4])
y1 <- y0 + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
segments(t0, y0, t0 + h, y1, lwd = 2, lty = 2, col = oi["green"])
points(t0 + h, y1, pch = 15, cex = 1.8, col = oi["green"])
points(t0 + h, exact(t0 + h), pch = 4, cex = 2, lwd = 3)
legend("topleft", c("k1 at the start", "k2 at the half step (via k1)",
                    "k3 at the half step (via k2)", "k4 at the end (via k3)",
                    "RK4 step: weights 1, 2, 2, 1 (/6)", "x  exact value"),
       col = c(cols, oi["green"], "black"), lwd = c(4, 4, 4, 4, 2, NA),
       lty = c(1, 1, 1, 1, 2, NA), cex = 1.2, bty = "n")
dev.off()

## ---- fig04 + fig05: order and cost ----------------------------------------

hs <- 2^(1:-8)
meth <- list(Euler = euler_step, Heun = heun_step, RK4 = rk4_step)
stages <- c(1, 2, 4)
err <- sapply(meth, function(s)
  sapply(hs, function(h) abs(rk_solve(s, flog, 0, N0, h, 10 / h)$y[10 / h + 1] - exact(10))))
mcols <- oi[c("red", "orange", "blue")]

svg("fig04-order.svg", width = 9, height = 5.8)
par(mar = c(5, 6, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(hs, err[, 1], log = "xy", type = "n", ylim = c(1e-14, 0.1), xlim = c(0.004, 2),
     xlab = "step size h (hours)", ylab = "", axes = FALSE)
title(ylab = "error in N(10)", line = 4.5)
axis(1, at = c(0.01, 0.1, 1), labels = c("0.01", "0.1", "1"))
axis(2, at = 10^c(-14, -11, -8, -5, -2), labels = c("1e-14", "1e-11", "1e-8", "1e-5", "0.01"), las = 1); box()
for (j in 1:3) lines(hs, err[, j], type = "o", pch = 19, lwd = 3, cex = 1.2, col = mcols[j])
for (j in 1:3) {
  p <- c(1, 2, 4)[j]; i0 <- 5
  lines(hs[4:10], err[i0, j] * (hs[4:10] / hs[i0])^p * 3, lty = 2, lwd = 2, col = oi["grey"])
}
text(0.006, c(4e-4, 1.2e-6, 2e-12), c("slope 1", "slope 2", "slope 4"), cex = 1.3, col = oi["grey"])
legend("bottomright", c("Euler", "Heun", "RK4"), col = mcols, lwd = 3, pch = 19, cex = 1.3, bty = "n")
dev.off()

svg("fig05-cost.svg", width = 9, height = 5.8)
par(mar = c(5, 6, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
nf <- outer(10 / hs, stages)
plot(nf[, 1], err[, 1], log = "xy", type = "n", ylim = c(1e-14, 0.1), xlim = c(5, 2e4),
     xlab = "number of f evaluations (t from 0 to 10)", ylab = "", axes = FALSE)
title(ylab = "error in N(10)", line = 4.5)
axis(1, at = c(10, 100, 1000, 10000), labels = c("10", "100", "1000", "10 000"))
axis(2, at = 10^c(-14, -11, -8, -5, -2), labels = c("1e-14", "1e-11", "1e-8", "1e-5", "0.01"), las = 1); box()
for (j in 1:3) lines(nf[, j], err[, j], type = "o", pch = 19, lwd = 3, cex = 1.2, col = mcols[j])
abline(h = 1e-6, lty = 3, lwd = 2, col = oi["grey"])
text(8, 3e-7, "1e-6", cex = 1.2, col = oi["grey"])
legend("bottomleft", c("Euler (1 per step)", "Heun (2 per step)", "RK4 (4 per step)"),
       col = mcols, lwd = 3, pch = 19, cex = 1.3, bty = "n")
dev.off()

## ---- fig06: repressilator and adaptive steps -------------------------------

s <- dp5(repress, 0, rep0, 50, rtol = 1e-6, atol = 1e-6)
svg("fig06-repressilator.svg", width = 12, height = 6)
layout(matrix(1:2, 2), heights = c(1.15, 1))
par(mar = c(1, 6, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
pc <- oi[c("blue", "red", "green")]
tq <- seq(0, 50, by = 0.05)
plot(range(s$t), c(0, max(s$y[, 4:6])), type = "n", xaxt = "n", las = 1,
     xlab = "", ylab = "protein p")
for (j in 4:6) lines(s$t, s$y[, j], lwd = 3, col = pc[j - 3])
legend("topright", c("LacI", "TetR", "cI"), col = pc, lwd = 3, cex = 1.2, bty = "n", horiz = TRUE)
par(mar = c(5, 6, 0.5, 1))
plot(s$t[-1], s$h, log = "y", type = "h", lwd = 1, col = oi["grey"], las = 1,
     xlab = "time (mRNA lifetimes)", ylab = "step size h")
points(s$t[-1], s$h, pch = 19, cex = 0.5, col = oi["orange"])
abline(h = 0.01, lty = 2, lwd = 2, col = oi["purple"])
text(44, 0.013, "fixed RK4 needs h = 0.01", col = oi["purple"], cex = 1.2, pos = 3)
dev.off()

## ---- fig07: conserved quantity of the predator-prey model --------------------

svg("fig07-lv-invariant.svg", width = 10, height = 5.6)
par(mar = c(5, 6, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
runs <- list(list(euler_step, 0.01, oi["red"], "Euler, h = 0.01"),
             list(heun_step, 0.01, oi["orange"], "Heun, h = 0.01"),
             list(rk4_step, 0.1, oi["blue"], "RK4, h = 0.1"),
             list(rk4_step, 0.01, oi["sky"], "RK4, h = 0.01"))
V0 <- V(lv0[1], lv0[2])
plot(NA, xlim = c(0, 50), ylim = c(1e-13, 1), log = "y", las = 1, axes = FALSE,
     xlab = "time (years)", ylab = "")
title(ylab = "|V(t) - V(0)|", line = 4.5)
axis(1); axis(2, at = 10^c(-12, -9, -6, -3, 0), labels = c("1e-12", "1e-9", "1e-6", "0.001", "1"), las = 1); box()
for (ru in runs) {
  so <- rk_solve(ru[[1]], lv, 0, lv0, ru[[2]], round(50 / ru[[2]]))
  dv <- abs(V(so$y[, 1], so$y[, 2]) - V0)
  keep <- seq(2, length(dv), by = max(1, round(0.1 / ru[[2]])))
  lines(so$t[keep], pmax(dv[keep], 1e-16), lwd = 2.5, col = ru[[3]])
}
legend("bottomright", sapply(runs, `[[`, 4), col = sapply(runs, `[[`, 3), lwd = 3,
       cex = 1.2, bty = "n")
dev.off()

## ---- fig08: tolerance: achieved error vs rtol, and the atol floor ------------

I0 <- 1e-6; s0 <- c(1 - I0, I0, 0)
ref <- dp5(sir, 0, s0, 60, rtol = 1e-13, atol = 1e-16)
Iref <- ref$y[nrow(ref$y), 2]
rts <- 10^-(2:10)
ach <- sapply(c(1e-6, 1e-16), function(at) sapply(rts, function(rt) {
  so <- dp5(sir, 0, s0, 60, rtol = rt, atol = at)
  abs(so$y[nrow(so$y), 2] - Iref) / Iref
}))
svg("fig08-tolerance.svg", width = 9, height = 5.8)
par(mar = c(5, 6, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(rts, ach[, 2], log = "xy", type = "o", pch = 19, lwd = 3, cex = 1.2, col = oi["blue"],
     ylim = c(1e-11, 1), xlab = "rtol", ylab = "", axes = FALSE)
title(ylab = "relative error in I(60)", line = 4.5)
axis(1, at = 10^c(-10, -8, -6, -4, -2), labels = c("1e-10", "1e-8", "1e-6", "1e-4", "0.01"))
axis(2, at = 10^c(-10, -8, -6, -4, -2, 0), labels = c("1e-10", "1e-8", "1e-6", "1e-4", "0.01", "1"), las = 1); box()
lines(rts, ach[, 1], type = "o", pch = 17, lwd = 3, cex = 1.3, col = oi["red"])
lines(rts, rts, lty = 2, lwd = 2, col = oi["grey"])
text(3e-9, 1.2e-9, "error = rtol", col = oi["grey"], cex = 1.2, pos = 4)
legend("topleft", c("atol = 1e-16", "atol = 1e-6  (= I at t = 0)"), col = oi[c("blue", "red")],
       pch = c(19, 17), lwd = 3, cex = 1.3, bty = "n")
dev.off()

## ---- fig09: dense output ---------------------------------------------------

y0s <- c(0.999, 0.001, 0)
s <- dp5(sir, 0, y0s, 120, rtol = 1e-4, atol = 1e-8)
svg("fig09-dense-output.svg", width = 11, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
fine <- dp5(sir, 0, y0s, 120, rtol = 1e-10, atol = 1e-12)
plot(fine$t, fine$y[, 2], type = "l", lwd = 6, col = oi["grey"], las = 1,
     xlab = "time (days)", ylab = "infected fraction I")
lines(s$t, s$y[, 2], lwd = 2, lty = 2, col = oi["red"])
tq <- seq(0, 120, by = 0.25)
iq <- findInterval(tq, s$t, rightmost.closed = TRUE)
dq <- hermite(tq, s$t[iq], s$t[iq + 1], s$y[iq, 2], s$y[iq + 1, 2], s$f[iq, 2], s$f[iq + 1, 2])
lines(tq, dq, lwd = 2.5, col = oi["blue"])
points(s$t, s$y[, 2], pch = 19, cex = 1.4, col = oi["red"])
legend("topright", c("exact", sprintf("solver steps (%d steps), joined by lines", length(s$h)),
                     "dense output between the steps"),
       col = oi[c("grey", "red", "blue")], lwd = c(6, 2, 2.5), lty = c(1, 2, 1),
       pch = c(NA, 19, NA), cex = 1.25, bty = "n")
dev.off()

## ---- fig10: events -------------------------------------------------------

svg("fig10-events.svg", width = 12, height = 5.3)
par(mfrow = c(1, 2), mar = c(5, 5.8, 2.5, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.4)
g <- beta * fine$y[, 1] - gamma
plot(fine$t, g, type = "l", lwd = 4, col = oi["grey"], xlim = c(20, 60), ylim = c(-0.1, 0.2),
     las = 1, xlab = "time (days)", ylab = expression(g == beta * S - gamma),
     main = "peak of I: g changes sign")
abline(h = 0, lty = 3)
gs <- beta * s$y[, 1] - gamma
points(s$t, gs, pch = 19, cex = 1.4, col = oi["red"])
i <- which(gs[-1] < 0 & gs[-length(gs)] > 0)
tq <- seq(s$t[i], s$t[i + 1], length.out = 50)
lines(tq, beta * hermite(tq, s$t[i], s$t[i + 1], s$y[i, 1], s$y[i + 1, 1], s$f[i, 1], s$f[i + 1, 1]) - gamma,
      lwd = 3, col = oi["blue"])
tpk <- uniroot(function(t) beta * hermite(t, s$t[i], s$t[i + 1], s$y[i, 1], s$y[i + 1, 1],
                                          s$f[i, 1], s$f[i + 1, 1]) - gamma, c(s$t[i], s$t[i + 1]))$root
points(tpk, 0, pch = 4, cex = 2.2, lwd = 3)
text(tpk, -0.03, sprintf("t = %.2f", tpk), cex = 1.3, pos = 4)
legend("topright", c("solver steps", "dense output in the step", "root of g"),
       col = c(oi[c("red", "blue")], "black"), pch = c(19, NA, 4), lwd = c(NA, 3, 3), cex = 1.15, bty = "n")
# missed event: prey dips below a threshold between two coarse steps
k <- rk_solve(rk4_step, lv, 0, lv0, 0.01, 1500)
plot(k$t, k$y[, 1], type = "l", lwd = 4, col = oi["grey"], xlim = c(0, 15), las = 1,
     xlab = "time (years)", ylab = "prey x", main = "threshold x = 3: missed by big steps")
abline(h = 3, lty = 2, lwd = 2, col = oi["purple"])
tc <- seq(0, 15, by = 2.5)
points(tc, k$y[match(round(tc, 2), round(k$t, 2)), 1], pch = 19, cex = 1.5, col = oi["red"])
legend("topright", c("solution", "values at steps of 2.5 yr", "threshold"),
       col = c(oi[c("grey", "red", "purple")]), lwd = c(4, NA, 2), lty = c(1, NA, 2),
       pch = c(NA, 19, NA), cex = 1.15, bty = "n")
dev.off()

cat("done\n")
