# Figures for the diffusion-implicit/ unit ("Implicit schemes and the method of lines").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd diffusion-implicit && Rscript make-diffusion-implicit-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito (colour-blind
# safe). Takes a few seconds. Makes:
#   fig01-embryo-setup.svg   the embryo between a cold and a warm laminar stream,
#                            and the 1D model along its long axis
#   fig02-amplification.svg  growth factor per step of each Fourier mode:
#                            FTCS, backward Euler, Crank-Nicolson, exact
#   fig03-cn-oscillation.svg "stop the flow": a 20/27 C step inside the embryo,
#                            Crank-Nicolson vs backward Euler with dt = 0.05 s
#   fig04-order.svg          verification on the chamber: error vs dt (BE, CN)
#                            and error vs dx, log-log
#   fig05-embryo-steady.svg  steady temperature inside the embryo with the
#                            streams running, and the 10-90 % width
#   fig06-embryo-time.svg    how fast the step forms after the streams are
#                            switched on, and the stop-flow decay
#
# All numbers are computed here with the same functions as on the slides.

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- solvers (as on the slides) ------------------------------------------------

thomas <- function(a, b, c, d) {
  n <- length(d)
  for (i in 2:n) {
    w <- a[i] / b[i - 1]
    b[i] <- b[i] - w * c[i - 1]
    d[i] <- d[i] - w * d[i - 1]
  }
  x <- numeric(n)
  x[n] <- d[n] / b[n]
  for (i in (n - 1):1) x[i] <- (d[i] - c[i] * x[i + 1]) / b[i]
  x
}

## ---- the embryo (1D along the anterior-posterior axis) -------------------------
## nodes x_0 .. x_n on [0, L], no flux through the two tips (ghost nodes),
## heat exchange with the stream outside at rate 1/tau (kex = 0: streams off)

D_T  <- 1.4e-7 * 1e12             # thermal diffusivity of water, um^2/s
L    <- 500                        # embryo length, um
R    <- 100                        # embryo radius, um
tau  <- R^2 / (2.405^2 * D_T)      # exchange time, s (slowest radial mode)
lam  <- sqrt(D_T * tau)            # = R / 2.405, um
ne   <- 200; dxe <- L / ne; xe <- (0:ne) * dxe; me <- ne + 1
Tout <- ifelse(xe < L / 2, 20, ifelse(xe > L / 2, 27, 23.5))
lo <- c(0, rep(1, me - 2), 2); up <- c(2, rep(1, me - 2), 0); di <- rep(-2, me)
Kmul <- function(u) di * u + lo * c(0, u[-me]) + up * c(u[-1], 0)
emb_step <- function(u, dt, theta, kex) {         # theta = 1: BE, 1/2: CN
  r <- D_T * dt / dxe^2; k <- kex * dt
  rhs <- u + (1 - theta) * (r * Kmul(u) + k * (Tout - u)) + theta * k * Tout
  thomas(-theta * r * lo, 1 - theta * r * di + theta * k, -theta * r * up, rhs)
}
emb_run <- function(u, dt, nsteps, theta, kex) {
  for (s in seq_len(nsteps)) u <- emb_step(u, dt, theta, kex)
  u
}
T_steady <- Tout - 3.5 * sign(xe - L / 2) * exp(-abs(xe - L / 2) / lam)
stopflow_exact <- function(t, K = 4000) {          # cosine series, no-flux ends
  s <- rep(23.5, me)
  for (k in seq(1, K, by = 2))
    s <- s - 14 / (k * pi) * sin(k * pi / 2) * cos(k * pi * xe / L) * exp(-D_T * (k * pi / L)^2 * t)
  s
}

## ---- fig01: set-up ---------------------------------------------------------------

svg("fig01-embryo-setup.svg", width = 12, height = 5.2)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(NA, xlim = c(0, 12), ylim = c(0, 5.2), axes = FALSE, xlab = "", ylab = "", asp = 1)
# channel seen from above: flow from top to bottom, cold stream left, warm right
rect(0.3, 0.4, 3.15, 5.0, col = adjustcolor(oi["blue"], 0.25), border = NA)
rect(3.15, 0.4, 6.0, 5.0, col = adjustcolor(oi["red"], 0.25), border = NA)
segments(c(0.3, 6.0), 0.4, c(0.3, 6.0), 5.0, lwd = 3)
arrows(c(1.1, 1.9, 4.4, 5.2), 4.8, c(1.1, 1.9, 4.4, 5.2), 3.9, lwd = 2.5, length = 0.12)
text(1.7, 1.0, "20 °C", cex = 1.6, col = oi["blue"], font = 2)
text(4.6, 1.0, "27 °C", cex = 1.6, col = oi["red"], font = 2)
th <- seq(0, 2 * pi, length.out = 200)
polygon(3.15 + 1.75 * cos(th), 2.7 + 0.7 * sin(th), col = "grey92", border = "black", lwd = 2.5)
text(3.15, 2.7, "embryo ~500 µm", cex = 1.3)
text(3.15, 0.1, "channel ~850 µm wide, streams at 5 cm/s", cex = 1.2)
# the 1D model
arrows(6.4, 2.7, 7.1, 2.7, lwd = 3, length = 0.15)
segments(7.4, 2.7, 11.6, 2.7, lwd = 6, col = "grey40")
points(seq(7.4, 11.6, length.out = 15), rep(2.7, 15), pch = 19, cex = 1.1)
text(9.5, 3.35, "grid along the long axis x", cex = 1.3)
text(7.4, 2.1, "x = 0", cex = 1.25); text(11.6, 2.1, "x = L", cex = 1.25)
text(9.5, 1.3, "no flux through the tips;\nexchange with the stream outside", cex = 1.2)
dev.off()

## ---- fig02: amplification factors ------------------------------------------------

z <- 10^seq(-2, 3, length.out = 400)   # z = 4 r sin^2(theta/2) = dt * |eigenvalue|
svg("fig02-amplification.svg", width = 10, height = 5.2)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(z, exp(-z), log = "x", type = "l", lwd = 4, col = "black", ylim = c(-1.6, 1.05),
     xlab = expression(Delta * t ~ "x decay rate of the mode   (" * 4 * r ~ sin^2 * "(" * theta * "/2))"),
     ylab = "growth factor per step g", las = 1)
abline(h = c(-1, 0, 1), col = oi["grey"], lty = c(2, 1, 2))
lines(z, 1 - z, lwd = 4, col = oi["orange"])
lines(z, 1 / (1 + z), lwd = 4, col = oi["blue"])
lines(z, (1 - z / 2) / (1 + z / 2), lwd = 4, col = oi["red"])
abline(v = 2, lty = 3, lwd = 2)
text(2.3, -1.45, "FTCS limit", pos = 4, cex = 1.2)
legend("bottomleft", c("exact  exp(-z)", "FTCS  1 - z", "backward Euler  1/(1+z)",
                       "Crank-Nicolson  (1-z/2)/(1+z/2)"),
       lwd = 4, col = c("black", oi["orange"], oi["blue"], oi["red"]), cex = 1.2, bty = "n")
dev.off()

## ---- fig03: CN oscillation with step data ---------------------------------------

u0 <- Tout                                  # the step, held by the streams until t = 0
svg("fig03-cn-oscillation.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
for (th in c(0.5, 1)) {
  plot(xe, stopflow_exact(0.1), type = "l", lwd = 5, col = "grey60", ylim = c(19.5, 27.5),
       xlab = "position along the embryo x (µm)", ylab = "temperature (°C)", las = 1,
       main = if (th == 0.5) "Crank-Nicolson, dt = 0.05 s" else "backward Euler, dt = 0.05 s")
  lines(xe, u0, lty = 3, lwd = 2)
  cols <- c(oi["orange"], oi["red"], oi["purple"])
  for (k in 1:3) lines(xe, emb_run(u0, 0.05, k, th, 0), lwd = 3, col = cols[k])
  legend("topleft", c("exact, t = 0.1 s", "t = 0.05", "t = 0.10", "t = 0.15"),
         lwd = c(5, 3, 3, 3), col = c("grey60", cols), cex = 1.1, bty = "n")
}
dev.off()

## ---- fig05: embryo steady state -----------------------------------------------------

Tss <- emb_run(rep(20, me), 0.005, 40, 1, 1 / tau)
f <- approxfun(Tss, xe)
x10 <- f(20.7); x90 <- f(26.3)
svg("fig05-embryo-steady.svg", width = 10, height = 5.2)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(xe, Tout, type = "s", lwd = 3, lty = 2, col = oi["grey"], ylim = c(19.5, 27.5),
     xlab = "position along the embryo x (µm)", ylab = "temperature (°C)", las = 1)
rect(x10, 19, x90, 28, col = adjustcolor(oi["yellow"], 0.35), border = NA)
lines(xe, Tss, lwd = 4, col = oi["red"])
points(xe[seq(1, me, by = 8)], T_steady[seq(1, me, by = 8)], pch = 1, cex = 1.2)
text((x10 + x90) / 2, 27.3, sprintf("10-90 %%: %.0f µm", x90 - x10), cex = 1.3)
legend("bottomright", c("stream outside", "backward Euler, t = 0.2 s",
                        "formula (no tips)"),
       lwd = c(3, 4, NA), lty = c(2, 1, NA), pch = c(NA, NA, 1),
       col = c(oi["grey"], oi["red"], "black"), cex = 1.2, bty = "n")
dev.off()

## ---- fig06: time course ---------------------------------------------------------------

dt <- 0.001; tt <- seq(0, 1, by = dt)
iA <- which.min(abs(xe - 125)); iP <- which.min(abs(xe - 375))
on <- off <- matrix(NA, length(tt), 2)
u <- rep(20, me); v <- Tout
on[1, ] <- u[c(iA, iP)]; off[1, ] <- v[c(iA, iP)]
for (k in 2:length(tt)) {
  u <- emb_step(u, dt, 0.5, 1 / tau); v <- emb_step(v, dt, 1, 0)
  on[k, ] <- u[c(iA, iP)]; off[k, ] <- v[c(iA, iP)]
}
svg("fig06-embryo-time.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
plot(tt * 1000, on[, 2], type = "l", lwd = 4, col = oi["red"], ylim = c(19.5, 27.5),
     xlim = c(0, 100), xlab = "time after switch-on (ms)", ylab = "temperature (°C)",
     las = 1, main = "streams switched on (embryo at 20 °C)")
lines(tt * 1000, on[, 1], lwd = 4, col = oi["blue"])
legend("right", c("posterior, x = 375 µm", "anterior, x = 125 µm"), lwd = 4,
       col = oi[c("red", "blue")], cex = 1.15, bty = "n")
plot(tt, off[, 2], type = "l", lwd = 4, col = oi["red"], ylim = c(19.5, 27.5),
     xlab = "time after the flow stops (s)", ylab = "temperature (°C)", las = 1,
     main = "flow stopped (insulated embryo)")
lines(tt, off[, 1], lwd = 4, col = oi["blue"])
dev.off()

## ---- fig04: order of accuracy on the chamber -----------------------------------------
## (filled in below)
