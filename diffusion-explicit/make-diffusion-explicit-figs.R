# Figures for the diffusion-explicit/ unit
# ("The diffusion (heat) equation by finite differences").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd diffusion-explicit && Rscript make-diffusion-explicit-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito (colour-blind safe).
# Running example: a flow-free chemotaxis channel of length L = 500 um between a
# source reservoir (c = 1 at x = 0) and a sink reservoir (c = 0 at x = L),
# small molecule D = 1000 um^2/s (1e-9 m^2/s), grid N = 50 (dx = 10 um), r = 0.4.
# Makes:
#   fig01-chamber.svg      schematic of the chamber and the profile that forms in it
#   fig02-space-time.svg   c(x, t) as a heat map, and the two slices that give the partial derivatives
#   fig03-grid.svg         the space-time grid, boundary/initial values and the FTCS stencil
#   fig04-profiles.svg     FTCS profiles against the exact (Fourier series) solution
#   fig05-scaling.svg      distance from the linear gradient vs t, and vs t D / L^2 (collapse)
#   fig06-convergence.svg  error at t = 30 s vs dx (r fixed), log-log, slope 2
#   fig07-unstable.svg     r = 0.5 vs r = 0.52: the zigzag instability
#   fig08-modes.svg        eigenvalues of the second-difference matrix and the slowest/fastest modes
#   fig09-noflux.svg       open chamber (sink) vs dead-end chamber (no-flux wall)

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999", black = "#000000")
um <- "µm"

## ---- the running example (same code as on the slides) ----------------------

D <- 1000; L <- 500; N <- 50; dx <- L / N
x <- seq(0, L, by = dx)
r <- 0.4; dt <- r * dx^2 / D
ftcs <- function(c, r, nsteps) {
  i <- 2:(length(c) - 1)
  for (k in seq_len(nsteps)) c[i] <- c[i] + r * (c[i - 1] - 2 * c[i] + c[i + 1])
  c
}
ftcs_closed <- function(c, r, nsteps) {          # no-flux wall at x = L
  n <- length(c); i <- 2:(n - 1)
  for (k in seq_len(nsteps)) {
    cn <- c[n] + 2 * r * (c[n - 1] - c[n])
    c[i] <- c[i] + r * (c[i - 1] - 2 * c[i] + c[i + 1])
    c[n] <- cn
  }
  c
}
exact <- function(x, t, L, D, nterm = 200) {
  s <- 1 - x / L
  for (n in 1:nterm) s <- s - 2 / (n * pi) * sin(n * pi * x / L) * exp(-(n * pi / L)^2 * D * t)
  s
}
c0 <- c(1, rep(0, N))

par_big <- function(...) par(cex.lab = 1.5, cex.axis = 1.35, cex.main = 1.5, ...)

## ---- fig01: the chamber -----------------------------------------------------

svg("fig01-chamber.svg", width = 11, height = 6.2)
layout(matrix(1:2, 2), heights = c(1, 1.15))
par_big(mar = c(0.5, 5, 0.5, 1))
plot(NA, xlim = c(-260, 760), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "")
rect(-250, 0.05, -10, 0.95, col = "#0072B233", border = oi["blue"], lwd = 2)
rect(510, 0.05, 750, 0.95, col = "grey95", border = oi["grey"], lwd = 2)
rect(-10, 0.38, 510, 0.62, col = NA, border = "black", lwd = 2)
xs <- seq(0, 500, length.out = 300)
for (j in seq_along(xs)[-1])
  rect(xs[j - 1], 0.385, xs[j], 0.615, border = NA,
       col = adjustcolor(oi["blue"], alpha.f = 0.25 * (1 - xs[j] / 500)))
set.seed(3)
points(runif(9, 40, 460), runif(9, 0.43, 0.57), pch = 21, bg = oi["orange"], cex = 2.2)
text(-130, 0.5, "source\nreservoir\nc = 1", cex = 1.4)
text(630, 0.5, "sink\nreservoir\nc = 0", cex = 1.4)
arrows(0, 0.25, 500, 0.25, code = 3, length = 0.1, lwd = 1.5)
text(250, 0.16, paste0("channel, L = 500 ", um, ", no flow"), cex = 1.4)
arrows(150, 0.78, 350, 0.78, length = 0.12, lwd = 2, col = oi["blue"])
text(250, 0.87, "diffusive flux", col = oi["blue"], cex = 1.35)
par(mar = c(4.8, 5, 0.5, 1))
plot(NA, xlim = c(-260, 760), ylim = c(0, 1.05), axes = FALSE,
     xlab = paste0("position x (", um, ")"), ylab = "")
axis(1, at = seq(0, 500, 100)); axis(2, pos = -20, las = 1)
text(-185, 0.5, "c / c(source)", cex = 1.4)
cols <- c(oi["sky"], oi["green"], oi["orange"], oi["red"])
ts <- c(5, 20, 64, 1000)
for (j in seq_along(ts)) lines(xs, exact(xs, ts[j], L, D), lwd = 3, col = cols[j],
                               lty = if (j == 4) 2 else 1)
legend(560, 1.05, c("t = 5 s", "20 s", "64 s", "steady (linear)"), col = cols,
       lwd = 3, lty = c(1, 1, 1, 2), bty = "n", cex = 1.25)
dev.off()

## ---- fig02: function of space and time, two slices ---------------------------

png("fig02-space-time.png", width = 13, height = 4.8, units = "in", res = 150)  # PNG: dense image
par_big(mfrow = c(1, 3), mar = c(4.8, 5, 3, 1))
tt <- seq(0, 120, length.out = 241); xx <- seq(0, L, length.out = 201)
Z <- outer(xx, tt, function(a, b) exact(a, pmax(b, 1e-3), L, D))
image(xx, tt, Z, col = hcl.colors(60, "Blues 3", rev = TRUE),
      xlab = paste0("x (", um, ")"), ylab = "t (s)", main = "c(x, t)")
abline(h = 30, col = oi["orange"], lwd = 3); abline(v = 100, col = oi["red"], lwd = 3)
plot(xx, exact(xx, 30, L, D), type = "l", lwd = 3, col = oi["orange"],
     xlab = paste0("x (", um, ")"), ylab = "c", main = "slice t = 30 s: c(x, 30)")
x1 <- 100; c1 <- exact(x1, 30, L, D); s1 <- (exact(x1 + 1e-3, 30, L, D) - exact(x1 - 1e-3, 30, L, D)) / 2e-3
segments(x1 - 80, c1 - 80 * s1, x1 + 80, c1 + 80 * s1, lty = 2, lwd = 2)
points(x1, c1, pch = 19, cex = 1.6)
text(x1 + 60, c1 + 0.12, expression(partialdiff * c / partialdiff * x), cex = 2.3, pos = 4)
plot(tt, exact(100, pmax(tt, 1e-3), L, D), type = "l", lwd = 3, col = oi["red"],
     xlab = "t (s)", ylab = "c", main = paste0("slice x = 100 ", um, ": c(100, t)"))
t1 <- 30; s2 <- (exact(100, t1 + 1e-3, L, D) - exact(100, t1 - 1e-3, L, D)) / 2e-3
segments(t1 - 20, c1 - 20 * s2, t1 + 20, c1 + 20 * s2, lty = 2, lwd = 2)
points(t1, c1, pch = 19, cex = 1.6)
text(t1 + 5, c1 - 0.12, expression(partialdiff * c / partialdiff * t), cex = 2.3, pos = 4)
dev.off()

## ---- fig03: the grid and the stencil ------------------------------------------

svg("fig03-grid.svg", width = 10, height = 5.6)
par_big(mar = c(4.5, 5.5, 1, 1))
nx <- 8; nt <- 5
plot(NA, xlim = c(-0.5, nx + 0.5), ylim = c(-0.5, nt + 0.3), axes = FALSE,
     xlab = expression("space   " * x[i] == i * Delta * x), ylab = expression("time   " * t[n] == n * Delta * t))
axis(1, at = 0:nx, labels = c(expression(x[0]), paste0("x", 1:(nx - 1)), expression(x[N])), tick = FALSE)
axis(2, at = 0:nt, labels = paste0("t", 0:nt), las = 1, tick = FALSE)
abline(h = 0:nt, col = "grey85"); abline(v = 0:nx, col = "grey85")
g <- expand.grid(i = 0:nx, n = 0:nt)
points(g$i, g$n, pch = 21, bg = "white", cex = 1.6)
points(rep(0, nt + 1), 0:nt, pch = 22, bg = oi["blue"], cex = 2.2)
points(rep(nx, nt + 1), 0:nt, pch = 22, bg = oi["grey"], cex = 2.2)
points(1:(nx - 1), rep(0, nx - 1), pch = 21, bg = oi["yellow"], cex = 2)
sx <- 4; sn <- 2
points(c(sx - 1, sx, sx + 1), rep(sn, 3), pch = 21, bg = oi["orange"], cex = 2.6)
points(sx, sn + 1, pch = 21, bg = oi["red"], cex = 2.6)
arrows(c(sx - 1, sx, sx + 1), sn + 0.08, sx, sn + 0.88, length = 0.1, lwd = 2, col = oi["red"])
text(c(sx - 1.15, sx + 0.12, sx + 1.15), sn - 0.35, c("r", "1 - 2r", "r"), cex = 1.4)
text(sx + 0.15, sn + 1.35, expression(c[i]^{n + 1}), cex = 1.5)
text(-0.5, nt + 0.25, "source c = 1", col = oi["blue"], cex = 1.2, pos = 4)
text(nx + 0.5, nt + 0.25, "sink c = 0", col = "grey30", cex = 1.2, pos = 2)
text(nx / 2, -0.45, "initial values: channel empty", cex = 1.2, col = "grey30")
dev.off()

## ---- fig04: FTCS vs exact -------------------------------------------------------

svg("fig04-profiles.svg", width = 9, height = 5.4)
par_big(mar = c(4.8, 5, 1, 1))
plot(NA, xlim = c(0, L), ylim = c(0, 1), xlab = paste0("x (", um, ")"), ylab = "c")
ts <- c(5, 15, 30, 64.4, 120)
cols <- c(oi["sky"], oi["blue"], oi["green"], oi["orange"], oi["red"])
xf <- seq(0, L, length.out = 400)
for (j in seq_along(ts)) {
  lines(xf, exact(xf, ts[j], L, D), lwd = 2.5, col = cols[j])
  cc <- ftcs(c0, r, round(ts[j] / dt))
  points(x[seq(1, 51, 2)], cc[seq(1, 51, 2)], pch = 19, col = cols[j], cex = 1.1)
}
lines(c(0, L), c(1, 0), lty = 2, lwd = 2)
legend("topright", c(paste0("t = ", ts, " s"), "steady state"), col = c(cols, "black"),
       lwd = c(rep(2.5, 5), 2), lty = c(rep(1, 5), 2), bty = "n", cex = 1.25)
text(10, 0.06, "dots: FTCS, lines: exact", pos = 4, cex = 1.25)
dev.off()

## ---- fig05: time to gradient, collapse ------------------------------------------

devtrace <- function(D, L, tmax) {
  dx <- L / N; dt <- r * dx^2 / D; x <- seq(0, L, by = dx); lin <- 1 - x / L
  ns <- round(tmax / dt); every <- max(1, ns %/% 300)
  c <- c0; out <- matrix(NA, 0, 2)
  for (k in seq(every, ns, by = every)) {
    c <- ftcs(c, r, every); out <- rbind(out, c(k * dt, max(abs(c - lin))))
  }
  out
}
cases <- list(c(1000, 500), c(1000, 1000), c(100, 500), c(100, 1000))
labs <- c(paste0("small molecule, L = 500 ", um), paste0("small molecule, L = 1000 ", um),
          paste0("protein, L = 500 ", um), paste0("protein, L = 1000 ", um))
cols <- c(oi["sky"], oi["blue"], oi["orange"], oi["red"])
tr <- lapply(cases, function(p) devtrace(p[1], p[2], 1.2 * p[2]^2 / p[1]))
svg("fig05-scaling.svg", width = 13, height = 5)
par_big(mfrow = c(1, 2), mar = c(4.8, 5, 2.5, 1))
plot(NA, xlim = c(5, 12000), ylim = c(1e-3, 1), log = "xy", xlab = "t (s)",
     ylab = "max |c - linear|", main = "distance from the linear gradient")
for (j in 1:4) lines(tr[[j]], lwd = 3, col = cols[j])
abline(h = 0.05, lty = 2)
legend("bottomleft", labs, col = cols, lwd = 3, bty = "n", cex = 1.05)
plot(NA, xlim = c(0, 1.2), ylim = c(1e-3, 1), log = "y", xlab = expression(t * D / L^2),
     ylab = "max |c - linear|", main = "same curves against t D / L²")
for (j in 4:1) lines(tr[[j]][, 1] * cases[[j]][1] / cases[[j]][2]^2, tr[[j]][, 2],
                     lwd = 7 - j, col = cols[j])
abline(h = 0.05, lty = 2); abline(v = 0.2578, lty = 3)
text(0.27, 0.3, "0.258", pos = 4, cex = 1.3)
dev.off()

## ---- fig06: convergence -----------------------------------------------------------

Ns <- c(20, 40, 80, 160, 320); err <- numeric(5)
for (j in 1:5) {
  dxj <- L / Ns[j]; dtj <- 0.4 * dxj^2 / D; xj <- seq(0, L, by = dxj)
  cj <- ftcs(c(1, rep(0, Ns[j])), 0.4, round(30 / dtj))
  err[j] <- max(abs(cj - exact(xj, 30, L, D, 2000)))
}
svg("fig06-convergence.svg", width = 8, height = 5.4)
par_big(mar = c(4.8, 5.5, 1, 1))
plot(L / Ns, err, log = "xy", pch = 19, cex = 1.6, col = oi["blue"],
     xlab = paste0("\u0394x (", um, "),  with r = 0.4 fixed"),
     ylab = "max error at t = 30 s", xlim = c(0.7, 30), ylim = c(2e-6, 2e-3))
lines(L / Ns, err, lwd = 2, col = oi["blue"])
lines(L / Ns, err[1] * (Ns[1] / Ns)^2, lty = 2, lwd = 2)
text(1.4, 2e-4, "dashed: slope 2", cex = 1.4, pos = 4)
text(L / Ns, err, paste(round(30 / (0.4 * (L / Ns)^2 / D)), "steps"), pos = 2, cex = 1.1, offset = 0.8)
dev.off()

## ---- fig07: instability -------------------------------------------------------------

svg("fig07-unstable.svg", width = 13, height = 5)
par_big(mfrow = c(1, 2), mar = c(4.8, 5, 2.5, 1))
xz <- seq(0, 200, length.out = 300)
plot(NA, xlim = c(0, 200), ylim = c(-0.6, 1.25), xlab = paste0("x (", um, ")"), ylab = "c",
     main = "r = 0.50 (dt = 0.050 s): smooth")
abline(h = c(0, 1), col = "grey80")
for (ns in c(40, 60, 80)) {
  lines(x, ftcs(c0, 0.5, ns), type = "o", pch = 19, lwd = 2, col = oi[c("sky", "blue", "green")][ns / 20 - 1])
}
legend("topright", c("40 steps", "60 steps", "80 steps"), col = oi[c("sky", "blue", "green")],
       lwd = 2, pch = 19, bty = "n", cex = 1.25)
plot(NA, xlim = c(0, 200), ylim = c(-0.6, 1.25), xlab = paste0("x (", um, ")"), ylab = "c",
     main = "r = 0.52 (dt = 0.052 s): zigzag grows")
abline(h = c(0, 1), col = "grey80")
for (ns in c(40, 60, 80)) {
  lines(x, ftcs(c0, 0.52, ns), type = "o", pch = 19, lwd = 2, col = oi[c("yellow", "orange", "red")][ns / 20 - 1])
}
legend("topright", c("40 steps", "60 steps", "80 steps"), col = oi[c("yellow", "orange", "red")],
       lwd = 2, pch = 19, bty = "n", cex = 1.25)
dev.off()

## ---- fig08: eigenvalues and modes ---------------------------------------------------

n <- N - 1
A <- diag(-2, n); A[cbind(1:(n - 1), 2:n)] <- 1; A[cbind(2:n, 1:(n - 1))] <- 1
e <- eigen(A, symmetric = TRUE)
lam <- D / dx^2 * e$values
svg("fig08-modes.svg", width = 13, height = 5)
par_big(mfrow = c(1, 2), mar = c(4.8, 5.5, 2.5, 1))
plot(1:n, lam, pch = 19, col = oi["blue"], xlab = "mode k", ylab = expression(lambda[k] ~ "(per s)"),
     main = "eigenvalues of (D / dx²) A", ylim = c(-42, 0))
abline(h = -2 / 0.05, lty = 2, col = oi["red"], lwd = 2)
text(1, -35.5, "Euler limit -2/dt for dt = 0.05 s (r = 1/2)", pos = 4, col = oi["red"], cex = 1.2)
text(1, -24, "slowest: -0.039 per s", pos = 4, cex = 1.2)
text(49, -30, "fastest: -39.96", pos = 2, cex = 1.2)
xi <- x[2:N]
v1 <- e$vectors[, 1]; v1 <- v1 / max(abs(v1)) * sign(v1[25])
vn <- e$vectors[, n]; vn <- vn / max(abs(vn)) * sign(vn[1])
plot(xi, v1, type = "l", lwd = 3, col = oi["blue"], ylim = c(-1.1, 1.7),
     xlab = paste0("x (", um, ")"), ylab = "mode shape", main = "slowest (k = 1) and fastest (k = 49) mode")
lines(xi, vn, type = "o", pch = 19, cex = 0.7, lwd = 1.5, col = oi["red"])
lines(xi, v1, lwd = 4, col = oi["blue"])
legend("top", c("k = 1: half a sine, decays in 25 s", "k = 49: zigzag, decays in 0.025 s"),
       col = oi[c("blue", "red")], lwd = 2.5, bty = "n", cex = 1.15, bg = "white")
dev.off()

## ---- fig09: open vs dead-end chamber --------------------------------------------------

svg("fig09-noflux.svg", width = 13, height = 5)
par_big(mfrow = c(1, 2), mar = c(4.8, 5, 2.5, 1))
ts <- c(10, 30, 64.4, 150, 330)
cols <- c(oi["sky"], oi["blue"], oi["green"], oi["orange"], oi["red"])
plot(NA, xlim = c(0, L), ylim = c(0, 1), xlab = paste0("x (", um, ")"), ylab = "c",
     main = "sink at x = L (c = 0): a gradient")
for (j in seq_along(ts)) lines(x, ftcs(c0, r, round(ts[j] / dt)), lwd = 3, col = cols[j])
legend("topright", paste0("t = ", ts, " s"), col = cols, lwd = 3, bty = "n", cex = 1.2)
plot(NA, xlim = c(0, L), ylim = c(0, 1), xlab = paste0("x (", um, ")"), ylab = "c",
     main = "wall at x = L (no flux): fills up")
for (j in seq_along(ts)) lines(x, ftcs_closed(c0, r, round(ts[j] / dt)), lwd = 3, col = cols[j])
abline(v = L, lwd = 6, col = "grey40")
dev.off()
