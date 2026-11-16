# Figures for the stiffness/ unit ("Stiffness and implicit methods").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd stiffness && Rscript make-stiffness-figs.R
#
# Base R graphics only (grDevices::svg) plus deSolve for reference solutions.
# Colours are Okabe-Ito (colour-blind safe). Makes:
#   fig01-enzyme.svg           mass-action enzyme kinetics on a log time axis: fast binding, slow catalysis
#   fig02-qssa.svg             full mass-action model vs the Michaelis-Menten (QSSA) model
#   fig03-stability.svg        stability regions of Euler, RK4 and backward Euler in the z = h*lambda plane
#   fig04-rk45-steps.svg       step sizes of an adaptive Dormand-Prince RK45 on the enzyme model
#   fig05-euler-vs-beuler.svg  explicit Euler (h = 1 ms, 3 ms) vs backward Euler (h = 0.1 s) on the enzyme
#   fig06-order.svg            error vs h: explicit Euler, backward Euler, BDF2 on the enzyme
#   fig07-robertson.svg        Robertson's problem from 1e-6 to 1e11 s (lsoda)
#   fig08-damping.svg          backward Euler damps a 24-hour oscillation

library(deSolve)

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- the running examples (same code as on the slides) ---------------------

kf <- 100; kr <- 100; kcat <- 1                  # per uM per s, per s, per s
E0 <- 1; S0 <- 10                                # uM
Km <- (kr + kcat) / kf; Vmax <- kcat * E0        # 1.01 uM, 1 uM/s
f <- function(y) {                               # y = (S, E, C, P)
  b <- kf * y[2] * y[1] - kr * y[3]
  c(-b, -b + kcat * y[3], b - kcat * y[3], kcat * y[3])
}
jac <- function(y) rbind(c(-kf * y[2], -kf * y[1], kr,         0),
                         c(-kf * y[2], -kf * y[1], kr + kcat,  0),
                         c( kf * y[2],  kf * y[1], -kr - kcat, 0),
                         c( 0,          0,          kcat,      0))
enzyme <- function(t, y, p) list(f(y))
y0 <- c(S = S0, E = E0, C = 0, P = 0)

beuler <- function(f, jac, y0, h, n, tol = 1e-10) {
  Y <- matrix(NA, n + 1, length(y0)); Y[1, ] <- y <- y0
  I <- diag(length(y0))
  for (i in 1:n) {
    yold <- y                                    # Newton starts at y_n
    for (k in 1:20) {
      delta <- solve(I - h * jac(y), -(y - yold - h * f(y)))
      y <- y + delta
      if (sqrt(sum(delta^2)) < tol) break
    }
    Y[i + 1, ] <- y
  }
  Y
}
bdf2 <- function(f, jac, y0, y1, h, n, tol = 1e-10) {  # needs two start values
  Y <- matrix(NA, n + 1, length(y0)); Y[1, ] <- y0; Y[2, ] <- y1
  I <- diag(length(y0))
  for (i in 2:n) {
    y <- Y[i, ]
    for (k in 1:20) {
      G <- y - 4 / 3 * Y[i, ] + 1 / 3 * Y[i - 1, ] - 2 / 3 * h * f(y)
      delta <- solve(I - 2 / 3 * h * jac(y), -G)
      y <- y + delta
      if (sqrt(sum(delta^2)) < tol) break
    }
    Y[i + 1, ] <- y
  }
  Y
}
euler <- function(f, y0, h, n) {
  Y <- matrix(NA, n + 1, length(y0)); Y[1, ] <- y <- y0
  for (i in 1:n) { y <- y + h * f(y); Y[i + 1, ] <- y }
  Y
}

## ---- fig01: fast binding, slow catalysis --------------------------------------

tt <- c(0, 10^seq(-5, log10(30), length.out = 400))
ref <- ode(y0, tt, enzyme, NULL, method = "lsoda", rtol = 1e-10, atol = 1e-12)
svg("fig01-enzyme.svg", width = 10, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(ref[-1, "time"], ref[-1, "S"], log = "x", type = "l", lwd = 4, col = oi["blue"],
     ylim = c(0, 10.5), xlab = "time t (s), log scale", ylab = expression("concentration (" * mu * "M)"),
     las = 1, xaxt = "n")
axis(1, at = 10^(-5:1), labels = c("10 µs", "100 µs", "1 ms", "10 ms", "0.1 s", "1 s", "10 s"))
lines(ref[-1, "time"], ref[-1, "P"], lwd = 4, col = oi["green"])
lines(ref[-1, "time"], 10 * ref[-1, "C"], lwd = 4, col = oi["red"])
rect(1e-5, -1, 5e-3, 11, col = adjustcolor(oi["orange"], 0.12), border = NA)
text(7e-5, 6.2, "fast: binding\n(~1 ms)", cex = 1.4, col = oi["orange"])
text(1.2, 6.2, "slow: catalysis\n(~10 s)", cex = 1.4, col = "grey30")
legend(5e-3, 7.5, c("substrate S", "complex C  (x10)", "product P"), lwd = 4,
       col = oi[c("blue", "red", "green")], cex = 1.3, bty = "n")
dev.off()

## ---- fig02: full model vs QSSA -------------------------------------------------

tl <- seq(0, 30, by = 0.05)
full <- ode(y0, tl, enzyme, NULL, method = "lsoda", rtol = 1e-10, atol = 1e-12)
mm <- function(t, y, p) list(-Vmax * y / (Km + y))
qssa <- ode(c(S = S0), tl, mm, NULL, rtol = 1e-10, atol = 1e-12)
svg("fig02-qssa.svg", width = 10, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(tl, full[, "P"], type = "l", lwd = 4, col = oi["green"], ylim = c(0, 10.5),
     xlab = "time t (s)", ylab = expression("product P (" * mu * "M)"), las = 1)
lines(tl, S0 - qssa[, "S"], lwd = 3, lty = 2, col = "black")
legend("bottomright", c("full mass action (4 ODEs, stiff)",
                        "Michaelis-Menten / QSSA (1 ODE, not stiff)"),
       lwd = c(4, 3), lty = c(1, 2), col = c(oi["green"], "black"), cex = 1.3, bty = "n")
dev.off()

## ---- fig03: stability regions -----------------------------------------------------

xs <- seq(-5, 3, length.out = 401); ys <- seq(-4, 4, length.out = 401)
z <- outer(xs, ys, function(a, b) complex(real = a, imaginary = b))
R_euler <- Mod(1 + z)
R_rk4 <- Mod(1 + z + z^2 / 2 + z^3 / 6 + z^4 / 24)
R_be <- Mod(1 / (1 - z))
svg("fig03-stability.svg", width = 13, height = 5)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1), cex.lab = 1.6, cex.axis = 1.4, cex.main = 1.8)
panel <- function(R, main, col, outside = FALSE) {
  plot(NA, xlim = range(xs), ylim = range(ys), asp = 1, xlab = expression("Re " * z),
       ylab = expression("Im " * z), main = main, las = 1)
  fill <- adjustcolor(col, 0.45)
  if (outside) rect(-10, -10, 10, 10, col = fill, border = NA)   # stable outside the curve
  for (cl in contourLines(xs, ys, R, levels = 1))
    polygon(cl$x, cl$y, col = if (outside) "white" else fill, border = col, lwd = 2.5)
  abline(h = 0, v = 0, col = oi["grey"])
  box()
}
panel(R_euler, "Euler: |1 + z| ≤ 1", oi["blue"])
text(-2.3, -0.5, "-2", cex = 1.5)
panel(R_rk4, "RK4", oi["green"])
text(-3.4, -0.5, "-2.79", cex = 1.5)
panel(R_be, "backward Euler: |1/(1 - z)| ≤ 1", oi["red"], outside = TRUE)
text(-2.5, 2.5, "whole left\nhalf-plane", cex = 1.5)
dev.off()

## ---- fig04: adaptive RK45 on the enzyme ------------------------------------------
## Hand-written Dormand-Prince 5(4) with a standard step-size controller.

dp_a <- list(c(), c(1/5), c(3/40, 9/40), c(44/45, -56/15, 32/9),
             c(19372/6561, -25360/2187, 64448/6561, -212/729),
             c(9017/3168, -355/33, 46732/5247, 49/176, -5103/18656),
             c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84))
dp_b <- c(35/384, 0, 500/1113, 125/192, -2187/6784, 11/84, 0)
dp_e <- dp_b - c(5179/57600, 0, 7571/16695, 393/640, -92097/339200, 187/2100, 1/40)
rk45 <- function(f, y, tend, rtol = 1e-6, atol = 1e-6, h = 1e-4) {
  t <- 0; ts <- 0; hs <- numeric(0); nrej <- 0
  while (t < tend) {
    h <- min(h, tend - t)
    K <- matrix(0, 7, length(y))
    for (s in 1:7) {
      ys <- y
      if (s > 1) for (j in 1:(s - 1)) ys <- ys + h * dp_a[[s]][j] * K[j, ]
      K[s, ] <- f(ys)
    }
    ynew <- y + h * colSums(dp_b * K)
    err <- sqrt(mean((h * colSums(dp_e * K) / (atol + rtol * pmax(abs(y), abs(ynew))))^2))
    if (err <= 1) { t <- t + h; y <- ynew; ts <- c(ts, t); hs <- c(hs, h) } else nrej <- nrej + 1
    h <- h * min(5, max(0.2, 0.9 * err^(-1/5)))
  }
  list(t = ts, h = hs, nrej = nrej, y = y)
}
run <- rk45(f, unname(y0), 30)
cat("RK45 (hand-written): accepted", length(run$h), "rejected", run$nrej, "\n")
lam_fast <- sapply(seq(0, 30, by = 0.1), function(tk) {
  yk <- ode(y0, c(0, tk + 1e-12), enzyme, NULL, method = "lsoda", rtol = 1e-10, atol = 1e-12)[2, -1]
  min(Re(eigen(jac(yk))$values))
})
svg("fig04-rk45-steps.svg", width = 10, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(run$t[-1], run$h, log = "y", type = "l", col = oi["blue"], lwd = 1.5,
     ylim = c(1e-4, 2), xlab = "time t (s)", ylab = "step size h (s)", las = 1, yaxt = "n")
axis(2, at = 10^(-4:0), labels = c("1e-4", "0.001", "0.01", "0.1", "1"), las = 1)
lines(seq(0, 30, by = 0.1), 3.3 / abs(lam_fast), lwd = 3, lty = 2, col = oi["red"])
abline(h = 0.5, lwd = 3, lty = 3, col = oi["green"])
text(15, 0.3, "step the slow curve would allow for accuracy (~0.5 s)", cex = 1.3, col = oi["green"])
text(15, 1.5e-4, paste0("RK45: ", length(run$h), " accepted + ", run$nrej,
                        " rejected steps"), cex = 1.3, col = oi["blue"])
legend(0, 0.12, c("RK45 step size", "stability limit 3.3 / |fast eigenvalue|"),
       lwd = c(2, 3), lty = c(1, 2), col = oi[c("blue", "red")], cex = 1.3, bty = "n")
dev.off()

## ---- fig05: explicit vs backward Euler on the enzyme --------------------------------

e1 <- euler(f, unname(y0), 1e-3, 3000)
e3 <- euler(f, unname(y0), 3e-3, 1000)
be <- beuler(f, jac, unname(y0), 0.1, 30)
svg("fig05-euler-vs-beuler.svg", width = 10, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(tl, full[, "C"], type = "l", lwd = 5, col = oi["grey"], xlim = c(0, 3), ylim = c(-0.5, 1.5),
     xlab = "time t (s)", ylab = expression("complex C (" * mu * "M)"), las = 1)
lines((0:3000) * 1e-3, e1[, 3], lwd = 2.5, lty = 2, col = oi["blue"])
ok <- is.finite(e3[, 3]) & abs(e3[, 3]) < 5
lines(((0:1000) * 3e-3)[ok], e3[ok, 3], lwd = 2, col = oi["red"])
points((0:30) * 0.1, be[, 3], pch = 19, cex = 1.3, col = oi["green"])
legend("topright", c("reference", "Euler, h = 1 ms (3000 steps)", "Euler, h = 3 ms: blows up",
                     "backward Euler, h = 0.1 s (30 steps)"),
       lwd = c(5, 2.5, 2, NA), lty = c(1, 2, 1, NA), pch = c(NA, NA, NA, 19), col = oi[c("grey", "blue", "red", "green")],
       cex = 1.3, bty = "n")
dev.off()

## ---- fig06: error vs h ----------------------------------------------------------

## Start after the binding transient (t = 0.1 s, from the reference), integrate 10 s:
## the solution is smooth there but the problem is still stiff.
refy <- function(t) unname(ode(y0, c(0, t), enzyme, NULL, method = "lsoda", rtol = 1e-12, atol = 1e-14)[2, -1])
ys <- refy(0.1); Pref <- refy(10.1)[4]
lam_s <- min(Re(eigen(jac(ys))$values))
cat("fast eigenvalue at t = 0.1 s:", lam_s, " Euler limit 2/|lambda| =", 2 / abs(lam_s), "\n")
hs_imp <- c(0.5, 0.2, 0.1, 0.05, 0.02, 0.01)
err_be <- sapply(hs_imp, function(h) abs(beuler(f, jac, ys, h, round(10 / h))[round(10 / h) + 1, 4] - Pref))
err_b2 <- sapply(hs_imp, function(h)                  # second start value from the reference
  abs(bdf2(f, jac, ys, refy(0.1 + h), h, round(10 / h))[round(10 / h) + 1, 4] - Pref))
hs_ex <- c(0.0015, 0.001, 5e-4, 2e-4, 1e-4)
err_ex <- sapply(hs_ex, function(h) abs(euler(f, ys, h, round(10 / h))[round(10 / h) + 1, 4] - Pref))
print(rbind(hs_imp, err_be, err_b2)); print(rbind(hs_ex, err_ex))
svg("fig06-order.svg", width = 10, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(hs_imp, err_be, log = "xy", type = "o", pch = 19, lwd = 3, cex = 1.3, col = oi["red"],
     xlim = c(1e-4, 0.5), ylim = c(1e-8, 1), xlab = "step size h (s)",
     ylab = "error in P after 10 s", las = 1, xaxt = "n", yaxt = "n")
axis(1, at = 10^(-4:0), labels = c("1e-4", "0.001", "0.01", "0.1", "1"))
axis(2, at = 10^(-8:0), labels = parse(text = paste0("10^", -8:0)), las = 1)
lines(hs_imp, err_b2, type = "o", pch = 17, lwd = 3, cex = 1.3, col = oi["green"])
lines(hs_ex, err_ex, type = "o", pch = 15, lwd = 3, cex = 1.3, col = oi["blue"])
abline(v = 2 / abs(lam_s), lty = 2, col = oi["blue"], lwd = 2)
text(2 / abs(lam_s), 1e-7, "  Euler stability limit", adj = 0, cex = 1.2, col = oi["blue"])
legend("topleft", c("explicit Euler (slope 1, only left of the limit)",
                        "backward Euler (slope 1)", "BDF2 (slope 2)"),
       lwd = 3, pch = c(15, 19, 17), col = oi[c("blue", "red", "green")], cex = 1.2, bty = "n")
dev.off()

## ---- fig07: Robertson's problem ---------------------------------------------------

rober <- function(t, y, p) {
  r1 <- 0.04 * y[1]; r2 <- 3e7 * y[2]^2; r3 <- 1e4 * y[2] * y[3]
  list(c(-r1 + r3, r1 - r2 - r3, r2))
}
tr <- c(0, 10^seq(-6, 11, by = 0.05))
rob <- ode(c(y1 = 1, y2 = 0, y3 = 0), tr, rober, NULL, method = "lsoda", rtol = 1e-8, atol = 1e-12)
svg("fig07-robertson.svg", width = 10, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(rob[-1, "time"], rob[-1, "y1"], log = "x", type = "l", lwd = 4, col = oi["blue"],
     ylim = c(0, 1.05), xlab = "time t, log scale", ylab = "concentration", las = 1, xaxt = "n")
axis(1, at = 10^seq(-6, 10, by = 2), labels = parse(text = paste0("10^", seq(-6, 10, by = 2))))
lines(rob[-1, "time"], 1e4 * rob[-1, "y2"], lwd = 4, col = oi["red"])
lines(rob[-1, "time"], rob[-1, "y3"], lwd = 4, col = oi["green"])
legend("left", c(expression(y[1]), expression(y[2] %*% 10^4), expression(y[3])), lwd = 4,
       col = oi[c("blue", "red", "green")], cex = 1.4, bty = "n")
dev.off()

## ---- fig08: backward Euler damps an oscillation -------------------------------------

w <- 2 * pi / 24
A <- matrix(c(0, 1, -w^2, 0), 2, byrow = TRUE)
osc_f <- function(y) as.vector(A %*% y)
osc_j <- function(y) A
th <- seq(0, 96, by = 0.1)
b1 <- beuler(osc_f, osc_j, c(1, 0), 1, 96)
b01 <- beuler(osc_f, osc_j, c(1, 0), 0.1, 960)
svg("fig08-damping.svg", width = 10, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(th, cos(w * th), type = "l", lwd = 4, col = oi["grey"], ylim = c(-1.1, 1.4),
     xlab = "time t (h)", ylab = "clock protein deviation x", las = 1, xaxt = "n")
axis(1, at = seq(0, 96, by = 24))
lines(th, b01[, 1], lwd = 2.5, col = oi["orange"])
lines(0:96, b1[, 1], type = "o", pch = 19, cex = 0.8, lwd = 2.5, col = oi["red"])
legend("top", c("exact cos(ωt)", "backward Euler, h = 0.1 h", "backward Euler, h = 1 h"),
       lwd = c(4, 2.5, 2.5), pch = c(NA, NA, 19), col = oi[c("grey", "orange", "red")],
       cex = 1.3, bty = "n", horiz = FALSE, ncol = 3, x.intersp = 0.5)
dev.off()
