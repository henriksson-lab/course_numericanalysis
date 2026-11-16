# Figures for the euler/ unit ("Initial value problems and Euler's method").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd euler && Rscript make-euler-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Makes:
#   fig01-direction-field.svg  direction field of the logistic ODE with solutions
#   fig02-euler-steps.svg      Euler steps (h = 1 h): local vs global error
#   fig03-euler-vs-exact.svg   Euler with h = 1, 0.5, 0.1 h against the exact curve
#   fig04-order.svg            local error (slope 2) and global error (slope 1) vs h
#   fig05-unstable.svg         drug clearance C' = -0.8 C with h = 0.5, 2, 3 h
#   fig06-logistic-cycle.svg   logistic growth with h = 3 h: a spurious boom-bust cycle
#   fig07-sir.svg              SIR epidemic; error in the peak vs h
#   fig08-oscillator.svg       Euler on a 24-hour oscillator spirals outwards

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- the running examples (same code as on the slides) ---------------------

euler <- function(f, t0, y0, h, n) {            # scalar Euler, as on the slides
  t <- t0 + (0:n) * h
  y <- numeric(n + 1); y[1] <- y0
  for (i in 1:n) y[i + 1] <- y[i] + h * f(t[i], y[i])
  data.frame(t = t, y = y)
}
euler_sys <- function(f, t0, y0, h, n) {        # vector Euler, as on the slides
  t <- t0 + (0:n) * h
  Y <- matrix(0, n + 1, length(y0)); Y[1, ] <- y0
  for (i in 1:n) Y[i + 1, ] <- Y[i, ] + h * f(t[i], Y[i, ])
  list(t = t, y = Y)
}

K <- 1.2; N0 <- 0.05; r <- 0.8                 # logistic OD curve, per hour
f <- function(t, N) r * N * (1 - N / K)
exact_from <- function(t, t0, n0) K / (1 + (K / n0 - 1) * exp(-r * (t - t0)))
exact <- function(t) exact_from(t, 0, N0)

## ---- fig01: direction field --------------------------------------------------

svg("fig01-direction-field.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(NA, xlim = c(0, 10), ylim = c(0, 1.6), xlab = "time t (h)",
     ylab = "OD  N(t)", las = 1)
tg <- seq(0.25, 9.75, by = 0.5); Ng <- seq(0.05, 1.55, by = 0.1)
len <- 0.22                                       # segment length in plot units
sx <- 10 / 1.6 * 5.5 / 9                          # aspect correction (x per y unit)
for (tt in tg) for (NN in Ng) {
  s <- f(tt, NN)
  dx <- len / sqrt(1 + (s * sx)^2) ; dy <- s * dx
  segments(tt - dx / 2, NN - dy / 2, tt + dx / 2, NN + dy / 2, col = oi["grey"], lwd = 1.5)
}
tt <- seq(0, 10, length.out = 300)
for (n0 in c(0.01, 0.05, 0.2, 1.6)) lines(tt, exact_from(tt, 0, n0), lwd = 3,
                                          col = if (n0 == 0.05) oi["blue"] else oi["sky"])
abline(h = c(0, K), lty = 2, col = oi["red"], lwd = 2)
text(9.9, K + 0.07, "N = K = 1.2", adj = 1, cex = 1.3, col = oi["red"])
points(0, N0, pch = 19, cex = 1.6, col = oi["blue"])
text(0.25, N0 + 0.08, expression(N[0] == 0.05), adj = 0, cex = 1.3, col = oi["blue"])
dev.off()

## ---- fig02: Euler steps, local vs global error ------------------------------

svg("fig02-euler-steps.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
e1 <- euler(f, 0, N0, 1, 4)
tt <- seq(0, 4, length.out = 300)
plot(tt, exact(tt), type = "l", lwd = 3, col = oi["blue"], ylim = c(0, 0.65),
     xlab = "time t (h)", ylab = "OD  N(t)", las = 1)
for (i in 2:4) {                                  # exact curve through each Euler point
  ts <- seq(e1$t[i], 4, length.out = 200)
  lines(ts, exact_from(ts, e1$t[i], e1$y[i]), lty = 2, lwd = 2, col = oi["grey"])
}
lines(e1$t, e1$y, lwd = 3, col = oi["orange"])
points(e1$t, e1$y, pch = 19, cex = 1.5, col = oi["orange"])
# errors at t = 3
loc <- exact_from(3, 2, e1$y[3])
arrows(3.08, e1$y[4], 3.08, loc, code = 3, length = 0.08, lwd = 2, col = oi["green"])
text(3.15, (e1$y[4] + loc) / 2, "local error", adj = 0, cex = 1.3, col = oi["green"])
arrows(2.92, e1$y[4], 2.92, exact(3), code = 3, length = 0.08, lwd = 2, col = oi["red"])
text(2.85, (e1$y[4] + exact(3)) / 2 + 0.03, "global error", adj = 1, cex = 1.3, col = oi["red"])
legend("topleft", c("exact solution", "Euler, h = 1 h", "exact solution through an Euler point"),
       col = oi[c("blue", "orange", "grey")], lty = c(1, 1, 2), lwd = c(3, 3, 2),
       pch = c(NA, 19, NA), cex = 1.25, bty = "n")
dev.off()

## ---- fig03: Euler with three step sizes --------------------------------------

svg("fig03-euler-vs-exact.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
tt <- seq(0, 10, length.out = 300)
plot(tt, exact(tt), type = "l", lwd = 4, col = oi["blue"], ylim = c(0, 1.3),
     xlab = "time t (h)", ylab = "OD  N(t)", las = 1)
hs <- c(1, 0.5, 0.1); cols <- oi[c("red", "orange", "green")]
for (j in 1:3) {
  e <- euler(f, 0, N0, hs[j], round(10 / hs[j]))
  lines(e$t, e$y, lwd = 2.5, col = cols[j], type = if (hs[j] >= 0.5) "o" else "l", pch = 19)
}
legend("bottomright", c("exact", "Euler h = 1", "Euler h = 0.5", "Euler h = 0.1"),
       col = c(oi["blue"], cols), lwd = c(4, 2.5, 2.5, 2.5), pch = c(NA, 19, 19, NA),
       cex = 1.3, bty = "n")
dev.off()

## ---- fig04: local and global error vs h --------------------------------------

svg("fig04-order.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
h <- 2^-(0:10)
glob <- sapply(h, function(h) abs(tail(euler(f, 0, N0, h, 5 / h)$y, 1) - exact(5)))
loc  <- sapply(h, function(h) abs(N0 + h * f(0, N0) - exact(h)))
plot(h, glob, log = "xy", type = "o", pch = 19, lwd = 3, cex = 1.3, col = oi["red"],
     ylim = range(c(glob, loc)), xlab = "step size h (hours)", ylab = "error",
     axes = FALSE)
lines(h, loc, type = "o", pch = 17, lwd = 3, cex = 1.3, col = oi["green"])
axis(1, at = 10^(-3:0), labels = c("0.001", "0.01", "0.1", "1"))
axis(2, at = 10^(-9:0), labels = parse(text = paste0("10^", -9:0)), las = 1); box()
text(0.004, 1.6e-3, "global error at t = 5 h:  slope 1", adj = 0, cex = 1.3, col = oi["red"])
text(0.004, 1.2e-7, "local error, one step:  slope 2", adj = 0, cex = 1.3, col = oi["green"])
dev.off()

## ---- fig05: drug clearance, three step sizes ----------------------------------

svg("fig05-unstable.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
k <- 0.8; tt <- seq(0, 24, length.out = 300)
plot(tt, 10 * exp(-k * tt), type = "l", lwd = 4, col = oi["blue"], ylim = c(-15, 20),
     xlab = "time t (h)", ylab = "C (mg/L)", las = 1)
abline(h = 0, col = oi["grey"])
hs <- c(0.5, 2, 3); cols <- oi[c("green", "orange", "red")]
for (j in 1:3) {
  n <- floor(24 / hs[j]); tj <- (0:n) * hs[j]
  lines(tj, 10 * (1 - k * hs[j])^(0:n), type = "o", pch = 19, lwd = 2.5, cex = 1.2, col = cols[j])
}
legend("topright", c("exact", "h = 0.5:  factor 0.6", "h = 2:  factor -0.6", "h = 3:  factor -1.4"),
       col = c(oi["blue"], cols), lwd = c(4, 2.5, 2.5, 2.5), pch = c(NA, 19, 19, 19),
       cex = 1.3, bty = "n")
dev.off()

## ---- fig06: logistic with h = 3: spurious cycle ----------------------------

svg("fig06-logistic-cycle.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
tt <- seq(0, 40, length.out = 400)
plot(tt, exact(tt), type = "l", lwd = 4, col = oi["blue"], ylim = c(0, 1.6),
     xlab = "time t (h)", ylab = "OD  N(t)", las = 1)
e <- euler(f, 0, N0, 1, 40);   lines(e$t, e$y, type = "o", pch = 19, lwd = 2, col = oi["green"])
e <- euler(f, 0, N0, 3, 13);   lines(e$t, e$y, type = "o", pch = 19, lwd = 2.5, cex = 1.3, col = oi["red"])
legend("bottomright", c("exact", "Euler h = 1 h  (rh = 0.8)", "Euler h = 3 h  (rh = 2.4)"),
       col = oi[c("blue", "green", "red")], lwd = c(4, 2, 2.5), pch = c(NA, 19, 19),
       cex = 1.3, bty = "n")
dev.off()

## ---- fig07: SIR epidemic and the error in the peak ----------------------------

beta <- 0.3; gamma <- 0.1; Npop <- 10000
sir <- function(t, y) {
  S <- y[1]; I <- y[2]
  c(-beta * S * I / Npop, beta * S * I / Npop - gamma * I, gamma * I)
}
y0 <- c(9990, 10, 0)
Imax <- sum(y0[1:2]) - Npop / (beta / gamma) * (1 + log(beta / gamma * y0[1] / Npop))

svg("fig07-sir.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
ref <- euler_sys(sir, 0, y0, 0.01, 16000)
plot(ref$t, ref$y[, 1], type = "l", lwd = 3, col = oi["blue"], ylim = c(0, 10000),
     xlab = "time (days)", ylab = "people", las = 1)
lines(ref$t, ref$y[, 2], lwd = 3, col = oi["red"])
lines(ref$t, ref$y[, 3], lwd = 3, col = oi["green"])
e <- euler_sys(sir, 0, y0, 4, 40)
for (j in 1:3) points(e$t, e$y[, j], pch = 19, cex = 1, col = oi[c("blue", "red", "green")][j])
legend("right", c("S", "I", "R", "dots: Euler h = 4 d"), col = c(oi[c("blue", "red", "green")], "black"),
       lwd = c(3, 3, 3, NA), pch = c(NA, NA, NA, 19), cex = 1.25, bty = "n")
hs <- 2^(2:-6)
perr <- sapply(hs, function(h) max(euler_sys(sir, 0, y0, h, 200 / h)$y[, 2]) - Imax)
plot(hs, perr, log = "xy", type = "o", pch = 19, lwd = 3, cex = 1.3, col = oi["red"],
     xlab = "step size h (days)", ylab = "error in peak I (people)", axes = FALSE)
axis(1, at = c(0.01, 0.1, 1), labels = c("0.01", "0.1", "1"))
axis(2, at = c(1, 10, 100), las = 1); box()
lines(hs, perr[1] * hs / hs[1], lty = 2, lwd = 2, col = oi["grey"])
text(0.03, 30, "slope 1", cex = 1.3, col = oi["grey"])
dev.off()

## ---- fig08: Euler on a 24-hour oscillator --------------------------------------

w <- 2 * pi / 24
osc <- function(t, y) c(y[2], -w^2 * y[1])

svg("fig08-oscillator.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
e <- euler_sys(osc, 0, c(1, 0), 0.1, 2400)
plot(e$t / 24, e$y[, 1], type = "l", lwd = 2.5, col = oi["red"], ylim = c(-2.5, 2.5),
     xlab = "time (days)", ylab = "x (deviation)", las = 1)
tt <- seq(0, 240, length.out = 2000)
lines(tt / 24, cos(w * tt), lwd = 2.5, col = oi["blue"])
legend("topleft", c("exact  cos(wt)", "Euler h = 0.1 h"), col = oi[c("blue", "red")],
       lwd = 2.5, cex = 1.25, bty = "n", horiz = TRUE)
e <- euler_sys(osc, 0, c(1, 0), 1, 48)
plot(e$y[, 1], e$y[, 2] / w, type = "o", pch = 19, cex = 0.9, lwd = 2, col = oi["red"],
     asp = 1, xlab = "x", ylab = expression(x * "'" / omega), las = 1)
th <- seq(0, 2 * pi, length.out = 200)
lines(cos(th), sin(th), lwd = 3, col = oi["blue"])
points(1, 0, pch = 19, cex = 1.6, col = "black")
legend("topleft", c("exact: circle", "Euler h = 1 h, 2 days"), col = oi[c("blue", "red")],
       lwd = c(3, 2), pch = c(NA, 19), cex = 1.2, bty = "n")
dev.off()
