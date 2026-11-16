# Figures for the advection-diffusion/ unit ("Flow plus diffusion").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd advection-diffusion && Rscript make-advection-diffusion-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito (colour-blind
# safe). Takes about half a minute. Units: um and s throughout.
#
# Running examples (all parameters are our own round numbers, typical for
# microfluidics; all "data" are synthetic):
#   * 1D test problem: an inlet channel of length L = 500 um, flow u = 100 um/s
#     towards a culture chamber that holds a drug (c = 1 at x = L), clean medium
#     upstream (c = 0 at x = 0), D = 1000 um^2/s. Pe = 50.
#   * a 1 s pulse of a protein stimulus (D = 100 um^2/s) carried at u = 200 um/s
#     to cells 2 mm downstream.
#   * T-sensor / H-filter: width w = 100 um, mean speed 1000 um/s, Poiseuille
#     profile across the width, two streams side by side.
#   * Christmas-tree gradient generator: 2 (or 3) inlets, n -> n + 1 streams per
#     stage, 8 outlet streams of 50 um in a 400 um channel, D = 100, mean 200 um/s.
# Makes:
#   fig01-tsensor.svg      T-sensor schematic and the computed concentration
#                          fields of a small molecule and a protein
#   fig02-layer.svg        back-diffusion against the flow: exact profiles
#   fig03-wiggles.svg      central differences, cell Peclet 5, 2, 1, vs exact
#   fig04-upwind.svg       upwind differences vs exact, and central with D + u dx/2
#   fig05-order.svg        max error vs dx, central and upwind, log-log
#   fig06-pulse.svg        arrival of the pulse at the cells: numerical diffusion
#   fig07-cfl.svg          explicit upwind with C + 2r = 0.99 vs 1.045
#   fig08-poiseuille.svg   the given flow profile and the "local time" x / u(y)
#   fig09-march.svg        marching check: error vs dx (plug flow), and the
#                          explicit-marching blow-up at the walls
#   fig10-hfilter.svg      fraction crossed vs distance: small molecule, protein, cells
#   fig11-tree.svg         Christmas-tree network (2 inlets) and its outlet steps
#   fig12-outlet.svg       outlet profiles downstream: linear and hump gradients
#   fig13-fit.svg          synthetic T-sensor profile, fitted model and SSE(D)

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999", black = "#000000")
um <- "µm"
par_big <- function(...) par(cex.lab = 1.5, cex.axis = 1.35, cex.main = 1.5, ...)

## ---- solvers (same code as on the slides) ---------------------------------------

thomas <- function(a, b, c, d) {                 # from "Implicit schemes ..."
  n <- length(d)
  for (i in 2:n) { w <- a[i] / b[i - 1]; b[i] <- b[i] - w * c[i - 1]; d[i] <- d[i] - w * d[i - 1] }
  x <- numeric(n); x[n] <- d[n] / b[n]
  for (i in (n - 1):1) x[i] <- (d[i] - c[i] * x[i + 1]) / b[i]
  x
}

## 1D steady test problem
L <- 500; D <- 1000; u <- 100
exact1 <- function(x, u = 100, D = 1000) expm1(u * x / D) / expm1(u * L / D)
steady <- function(N, scheme, u = 100, D = 1000) {
  dx <- L / N; P <- u * dx / D; m <- N - 1          # P = cell Peclet number
  if (scheme == "central") { lo <- 1 + P / 2; di <- -2;     up <- 1 - P / 2 }
  else                     { lo <- 1 + P;     di <- -2 - P; up <- 1 }
  A <- diag(di, m)
  A[cbind(2:m, 1:(m - 1))] <- lo; A[cbind(1:(m - 1), 2:m)] <- up
  b <- numeric(m); b[m] <- -up                       # c(L) = 1
  c(0, solve(A, b), 1)
}

## pulse: explicit upwind in a long channel
up_pulse <- function(C, tend, D = 100, u = 200, dx = 10, xobs = 2000, Lc = 3000, tp = 1) {
  x <- seq(0, Lc, by = dx); n <- length(x); i <- 2:(n - 1)
  dt <- C * dx / u; r <- D * dt / dx^2
  c <- numeric(n); ns <- round(tend / dt); obs <- numeric(ns); iobs <- which.min(abs(x - xobs))
  snaps <- list()
  for (k in 1:ns) {
    old <- c
    c[i] <- old[i] - C * (old[i] - old[i - 1]) + r * (old[i - 1] - 2 * old[i] + old[i + 1])
    c[n] <- c[n - 1]                                  # outflow
    c[1] <- if (k * dt <= tp) 1 else 0                # inlet: 1 s pulse
    obs[k] <- c[iobs]
  }
  list(t = (1:ns) * dt, obs = obs, x = x, c = c, r = r)
}
erfc <- function(z) 2 * pnorm(-z * sqrt(2))
step_exact <- function(x, t, D, u = 200) {           # semi-infinite channel, step inflow
  if (t <= 0) return(0 * x)
  s <- 2 * sqrt(D * t)
  0.5 * erfc((x - u * t) / s) +
    0.5 * exp(u * x / D + log(2) + pnorm(-(x + u * t) / s * sqrt(2), log.p = TRUE))
}
pulse_exact <- function(x, t, D, tp = 1) step_exact(x, t, D) - step_exact(x, t - tp, D)

## T-sensor: marching in x, backward Euler, cell-centred grid across the width
w <- 100; ubar <- 1000
pois <- function(y, w = 100, ubar = 1000) 6 * ubar * (y / w) * (1 - y / w)
march <- function(c, u, D, dy, dx, nsteps) {
  M <- length(c); a <- D * dx / (u * dy^2)
  di <- 1 + 2 * a; di[c(1, M)] <- 1 + a[c(1, M)]     # no-flux walls
  for (k in seq_len(nsteps)) c <- thomas(-a, di, -a, c)
  c
}
march_explicit <- function(c, u, D, dy, dx, nsteps) {
  M <- length(c); a <- D * dx / (u * dy^2)
  for (k in seq_len(nsteps)) {
    g <- c(c[1], c, c[M])                             # mirror ghost cells
    c <- c + a * (g[1:M] - 2 * c + g[3:(M + 2)])
  }
  c
}
plug_exact <- function(y, x, D, K = 2001) {
  s <- 0.5
  for (n in 1:K) s <- s + 2 / (n * pi) * sin(n * pi / 2) * cos(n * pi * y / w) *
      exp(-(n * pi)^2 * D * x / (ubar * w^2))
  s
}
M <- 50; dy <- w / M; y <- (1:M - 0.5) * dy; uy <- pois(y)
c0 <- ifelse(y < w / 2, 1, 0)

## Christmas tree
tree_stage <- function(c) {                          # n streams -> n + 1 streams
  n <- length(c); k <- 0:n
  (k / n) * c(0, c) + ((n - k) / n) * c(c, 0)
}
tree <- function(inlets, nout) { c <- inlets; while (length(c) < nout) c <- tree_stage(c); c }

## ---- fig01: T-sensor ----------------------------------------------------------------

field <- function(D, xmax, dxs = 50) {
  cols <- matrix(NA, M, xmax / dxs + 1); c <- c0; cols[, 1] <- c
  for (k in 1:(xmax / dxs)) { c <- march(c, uy, D, dy, dxs, 1); cols[, k + 1] <- c }
  cols
}
F1 <- field(1000, 10000); F2 <- field(100, 10000)
xs <- seq(0, 10, by = 0.05)
svg("fig01-tsensor.svg", width = 13, height = 6.2)
layout(matrix(c(1, 1, 2, 3), 2), widths = c(1, 1.5))
par_big(mar = c(1, 1, 1, 1))
plot(NA, xlim = c(0, 10), ylim = c(0, 10), axes = FALSE, xlab = "", ylab = "", asp = 1)
polygon(c(0.3, 3.4, 3.4, 0.3), c(8.6, 6.4, 5.4, 7.6), col = "grey95", border = "black", lwd = 2)
polygon(c(0.3, 3.4, 3.4, 0.3), c(2.2, 4.4, 5.4, 3.2), col = "#0072B240", border = "black", lwd = 2)
polygon(c(3.4, 9.7, 9.7, 3.4), c(4.4, 4.4, 6.2, 5.4), col = "#0072B230", border = NA)
rect(3.4, 4.4, 9.7, 5.4, col = "#0072B230", border = NA)
rect(3.4, 4.4, 9.7, 6.4, col = NA, lwd = 2)
segments(3.4, 5.4, 9.7, 5.4, lty = 2, col = "grey40")
text(1.0, 8.8, "buffer", cex = 1.3, pos = 4)
text(0.6, 1.3, "sample: small molecules,\nproteins, cells", cex = 1.15, pos = 4)
arrows(4, 3.6, 9, 3.6, length = 0.12, lwd = 2); text(6.5, 3.0, "flow, x", cex = 1.3)
arrows(9.9, 4.4, 9.9, 6.4, code = 3, length = 0.08); text(9.9, 6.9, "w", cex = 1.3)
set.seed(2); points(runif(8, 3.7, 9.4), runif(8, 4.55, 5.25), pch = 21, bg = oi["orange"], cex = 1.8)
text(6.5, 7.4, "small molecules\nspread across", cex = 1.1, col = oi["blue"])
text(5, 9.4, "T-sensor / H-filter, top view", cex = 1.3, font = 2)
pal <- hcl.colors(60, "Blues 3", rev = TRUE)
par(mar = c(4.5, 5, 2.5, 1))
image(xs, y, t(F1), col = pal, zlim = c(0, 1), xlab = "", ylab = paste0("y (", um, ")"),
      main = paste0("small molecule, D = 1000 ", um, "²/s"))
image(xs, y, t(F2), col = pal, zlim = c(0, 1), xlab = "distance downstream x (mm)",
      ylab = paste0("y (", um, ")"), main = paste0("protein, D = 100 ", um, "²/s"))
dev.off()

## ---- fig02: back-diffusion boundary layer -----------------------------------------

svg("fig02-layer.svg", width = 12, height = 4.8)
layout(matrix(1:2, 1), widths = c(1, 1.25))
par_big(mar = c(1, 1, 1, 1))
plot(NA, xlim = c(0, 10), ylim = c(0, 6), axes = FALSE, xlab = "", ylab = "")
rect(0.3, 2.5, 6.5, 3.5, lwd = 2, col = "grey97")
rect(6.5, 0.8, 9.7, 5.2, lwd = 2, col = "#0072B240")
arrows(1.2, 3, 4.5, 3, lwd = 2.5, length = 0.12); text(2.8, 3.9, "flow u", cex = 1.4)
arrows(6.2, 2.2, 4.4, 2.2, lwd = 2, col = oi["blue"], length = 0.1)
text(5.3, 1.6, "diffusion", col = oi["blue"], cex = 1.3)
text(8.1, 3, "chamber\nwith drug\nc = 1", cex = 1.3)
text(0.3, 4.6, "medium in\nc = 0", cex = 1.3, pos = 4)
text(3.4, 0.4, paste0("inlet channel, L = 500 ", um), cex = 1.3)
par(mar = c(4.8, 5, 1, 1))
xx <- seq(0, L, length.out = 1000)
plot(NA, xlim = c(300, 500), ylim = c(0, 1), xlab = paste0("x (", um, ")"), ylab = "c")
us <- c(1, 10, 100); cl <- c(oi["sky"], oi["green"], oi["red"])
for (j in 1:3) lines(xx, exact1(xx, u = us[j]), lwd = 3, col = cl[j])
legend("topleft", paste0("u = ", us, " ", um, "/s,  Pe = ", us * L / D), col = cl, lwd = 3,
       bty = "n", cex = 1.25)
dev.off()

## ---- fig03: wiggles --------------------------------------------------------------

svg("fig03-wiggles.svg", width = 13, height = 4.6)
par_big(mfrow = c(1, 3), mar = c(4.8, 5, 3, 1))
for (N in c(10, 25, 50)) {
  xg <- seq(0, L, length.out = N + 1); cc <- steady(N, "central")
  plot(xx, exact1(xx), type = "l", lwd = 2.5, xlim = c(250, 500), ylim = c(-0.5, 1),
       xlab = paste0("x (", um, ")"), ylab = "c",
       main = sprintf("central, Δx = %g, cell Pe = %g", L / N, u * (L / N) / D))
  abline(h = 0, col = "grey70")
  lines(xg, cc, type = "o", pch = 19, col = oi["red"], lwd = 2, cex = 1.4)
}
dev.off()

## ---- fig04: upwind -----------------------------------------------------------------

svg("fig04-upwind.svg", width = 13, height = 4.6)
par_big(mfrow = c(1, 3), mar = c(4.8, 5, 3, 1))
for (N in c(10, 25, 50)) {
  xg <- seq(0, L, length.out = N + 1); cu <- steady(N, "upwind")
  plot(xx, exact1(xx), type = "l", lwd = 2.5, xlim = c(250, 500), ylim = c(-0.5, 1),
       xlab = paste0("x (", um, ")"), ylab = "c",
       main = sprintf("upwind, Δx = %g, cell Pe = %g", L / N, u * (L / N) / D))
  abline(h = 0, col = "grey70")
  lines(xx, exact1(xx, D = D + u * (L / N) / 2), lty = 2, lwd = 2, col = oi["grey"])
  lines(xg, cu, type = "o", pch = 19, col = oi["blue"], lwd = 2, cex = 1.4)
  if (N == 10) legend("topleft", c("exact, D", "exact, D + uΔx/2", "upwind"),
                      col = c("black", oi["grey"], oi["blue"]), lty = c(1, 2, 1),
                      lwd = 2, pch = c(NA, NA, 19), bty = "n", cex = 1.2)
}
dev.off()

## ---- fig05: order -------------------------------------------------------------------

Ns <- c(10, 20, 25, 40, 50, 100, 200, 400, 800)
ec <- sapply(Ns, function(N) max(abs(steady(N, "central") - exact1(seq(0, L, length.out = N + 1)))))
eu <- sapply(Ns, function(N) max(abs(steady(N, "upwind") - exact1(seq(0, L, length.out = N + 1)))))
cat("fig05 errors\n"); print(signif(cbind(dx = L / Ns, Pe_cell = u * L / Ns / D, central = ec, upwind = eu), 3))
svg("fig05-order.svg", width = 8.5, height = 5.4)
par_big(mar = c(4.8, 5.5, 1, 1))
plot(L / Ns, ec, log = "xy", type = "o", pch = 19, cex = 1.5, lwd = 2, col = oi["red"],
     xlab = paste0("Δx (", um, ")"), ylab = "max error", ylim = c(1e-5, 1))
lines(L / Ns, eu, type = "o", pch = 17, cex = 1.5, lwd = 2, col = oi["blue"])
abline(v = 2 * D / u, lty = 3); text(2 * D / u, 1.5e-5, "cell Pe = 2", pos = 2, cex = 1.2)
legend("topleft", c("central", "upwind"), col = c(oi["red"], oi["blue"]), pch = c(19, 17),
       lwd = 2, bty = "n", cex = 1.3)
dev.off()

## ---- fig06: pulse --------------------------------------------------------------------

tt <- seq(6, 16, by = 0.01)
svg("fig06-pulse.svg", width = 10, height = 5.4)
par_big(mar = c(4.8, 5, 1, 1))
plot(tt, sapply(tt, function(t) pulse_exact(2000, t, 100)), type = "l", lwd = 3,
     xlab = "time (s)", ylab = paste0("c at the cells (x = 2000 ", um, ")"), ylim = c(0, 1))
cl <- c(oi["blue"], oi["green"], oi["orange"])
for (j in 1:3) {
  dxj <- c(10, 5, 1)[j]; p <- up_pulse(0.5, 16, dx = dxj)
  lines(p$t, p$obs, lwd = 2.5, col = cl[j])
  cat("fig06 dx", dxj, "peak", max(p$obs), "\n")
}
lines(tt, sapply(tt, function(t) pulse_exact(2000, t, 600)), lty = 2, lwd = 2.5, col = oi["red"])
legend("topright", c(paste0("exact, D = 100 ", um, "²/s"),
                     paste0("upwind, Δx = ", c(10, 5, 1), " ", um, ", C = 0.5"),
                     "exact, D = 100 + 500"),
       col = c("black", cl, oi["red"]), lwd = c(3, 2.5, 2.5, 2.5, 2.5), lty = c(1, 1, 1, 1, 2),
       bty = "n", cex = 1.15)
dev.off()

## ---- fig07: CFL -------------------------------------------------------------------------

svg("fig07-cfl.svg", width = 13, height = 4.8)
par_big(mfrow = c(1, 2), mar = c(4.8, 5, 3, 1))
p1 <- up_pulse(0.9, 6); p2 <- up_pulse(0.95, 6)
cat("fig07 C+2r", 0.9 + 2 * p1$r, 0.95 + 2 * p2$r, "max", max(abs(p1$c)), max(abs(p2$c)), "\n")
plot(p1$x, p1$c, type = "l", lwd = 2.5, col = oi["blue"], xlim = c(0, 2000),
     xlab = paste0("x (", um, ")"), ylab = "c at t = 6 s",
     main = sprintf("C = 0.90, r = %.4f: C + 2r = %.3f", p1$r, 0.9 + 2 * p1$r))
plot(p2$x, p2$c, type = "l", lwd = 1.5, col = oi["red"], xlim = c(0, 2000),
     xlab = paste0("x (", um, ")"), ylab = "c at t = 6 s",
     main = sprintf("C = 0.95, r = %.4f: C + 2r = %.3f", p2$r, 0.95 + 2 * p2$r))
dev.off()

## ---- fig08: Poiseuille --------------------------------------------------------------------

svg("fig08-poiseuille.svg", width = 12, height = 4.8)
par_big(mfrow = c(1, 2), mar = c(4.8, 5.5, 1, 1))
yy <- seq(0, w, length.out = 400)
plot(pois(yy), yy, type = "l", lwd = 3, col = oi["blue"], xlab = paste0("u(y) (", um, "/s)"),
     ylab = paste0("y (", um, ")"), xlim = c(0, 1600))
abline(v = ubar, lty = 2); text(ubar, 8, "mean", pos = 4, cex = 1.2)
arrows(0, seq(10, 90, 10), pois(seq(10, 90, 10)), seq(10, 90, 10), length = 0.08, col = oi["grey"])
plot(1000 / pois(yy), yy, type = "l", lwd = 3, col = oi["red"], log = "x",
     xlab = paste0("time to reach x = 1 mm (s)"), ylab = paste0("y (", um, ")"), xlim = c(0.5, 100))
dev.off()

## ---- fig09: marching check + explicit failure -----------------------------------------

dxs <- c(200, 100, 50, 25, 12.5, 6.25)
em <- sapply(dxs, function(h) {
  Mm <- 200; dym <- w / Mm; ym <- (1:Mm - 0.5) * dym
  cm <- march(ifelse(ym < w / 2, 1, 0), rep(ubar, Mm), 1000, dym, h, 2000 / h)
  max(abs(cm - plug_exact(ym, 2000, 1000)))
})
cat("fig09 errors\n"); print(signif(cbind(dxs, em), 3))
growth <- function(dxm, xmax = 100) {               # max |c| along x, explicit marching
  c <- c0; out <- numeric(round(xmax / dxm))
  for (k in seq_along(out)) { c <- march_explicit(c, uy, 1000, dy, dxm, 1); out[k] <- max(abs(c)) }
  out
}
g1 <- growth(1); g01 <- growth(0.1)
cat("fig09 explicit limit near wall", min(uy) * dy^2 / (2 * 1000), "centre", max(uy) * dy^2 / (2 * 1000),
    " max|c| at x = 100: dx = 1:", g1[100], " dx = 0.1:", g01[1000], "\n")
svg("fig09-march.svg", width = 13, height = 4.8)
par_big(mfrow = c(1, 2), mar = c(4.8, 5.5, 3, 1))
plot(dxs, em, log = "xy", type = "o", pch = 19, cex = 1.5, lwd = 2, col = oi["blue"],
     xlab = paste0("\u0394x (", um, ")"), ylab = "max error at x = 2 mm",
     main = "implicit marching, plug flow")
lines(dxs, em[1] * dxs / dxs[1], lty = 2); text(30, em[1] * 30 / 200 * 0.6, "slope 1", pos = 1, cex = 1.3)
plot((1:100) * 1, g1, log = "y", type = "l", lwd = 3, col = oi["red"], ylim = c(0.5, 1e40),
     xlab = paste0("x (", um, ")"), ylab = "max |c| across the channel",
     main = "explicit marching, Poiseuille flow")
lines((1:1000) * 0.1, g01, lwd = 3, col = oi["blue"])
legend("topleft", paste0("\u0394x = ", c(1, 0.1), " ", um), col = c(oi["red"], oi["blue"]), lwd = 3,
       bty = "n", cex = 1.25)
dev.off()

## ---- fig10: H-filter -------------------------------------------------------------------------

right <- y > w / 2
xsteps <- seq(0, 20000, by = 100)
frac <- sapply(c(1000, 100, 0.05), function(Dm) {
  c <- c0; out <- 0
  for (k in 2:length(xsteps)) { c <- march(c, uy, Dm, dy, 100, 1); out <- c(out, sum(uy[right] * c[right]) / sum(uy * c0)) }
  out
})
cat("fig10 fraction crossed at 1, 2, 5, 10, 20 mm\n"); print(round(frac[xsteps %in% (c(1, 2, 5, 10, 20) * 1000), ], 4))
svg("fig10-hfilter.svg", width = 9, height = 5.4)
par_big(mar = c(4.8, 5.5, 1, 1))
plot(NA, xlim = c(0, 20), ylim = c(0, 0.52), xlab = "distance downstream (mm)",
     ylab = "fraction in the buffer half")
cl <- c(oi["blue"], oi["orange"], oi["red"])
for (j in 1:3) lines(xsteps / 1000, frac[, j], lwd = 3, col = cl[j])
abline(h = 0.5, lty = 2)
legend("right", c(paste0("small molecule, D = 1000 ", um, "²/s"),
                  paste0("protein, D = 100 ", um, "²/s"), paste0("cell, D ≈ 0.05 ", um, "²/s")),
       col = cl, lwd = 3, bty = "n", cex = 1.2)
dev.off()

## ---- fig11: Christmas tree --------------------------------------------------------------------

svg("fig11-tree.svg", width = 13, height = 5.6)
layout(matrix(1:2, 1), widths = c(1.35, 1))
par_big(mar = c(1, 1, 2.5, 1))
plot(NA, xlim = c(-0.5, 7.5), ylim = c(-0.8, 7.6), axes = FALSE, xlab = "", ylab = "",
     main = "2 inlets, 6 mixing stages, 8 outlets")
levels <- list(c(0, 1)); for (s in 1:6) levels[[s + 1]] <- tree_stage(levels[[s]])
pos <- function(n) (0:(n - 1)) - (n - 1) / 2 + 3.5
for (s in seq_along(levels)) {
  cs <- levels[[s]]; n <- length(cs); yl <- 7 - (s - 1) * 1.15
  if (s > 1) {
    pp <- pos(n - 1)
    for (k in 1:n) {
      if (k > 1) segments(pp[k - 1], yl + 1.15, pos(n)[k], yl, col = "grey60", lwd = 1.5)
      if (k < n) segments(pp[k], yl + 1.15, pos(n)[k], yl, col = "grey60", lwd = 1.5)
    }
  }
  points(pos(n), rep(yl, n), pch = 22, cex = 3.6, bg = gray(1 - 0.8 * cs))
  if (s == 7) text(pos(n), yl - 0.75, sprintf("%.2f", cs), cex = 1.05)
}
par(mar = c(4.8, 5, 2.5, 1))
yo <- seq(0, 400, length.out = 801)
step_prof <- function(cs, yv) cs[pmin(length(cs), floor(yv / (400 / length(cs))) + 1)]
plot(yo, step_prof(levels[[7]], yo), type = "s", lwd = 3, col = oi["blue"], ylim = c(0, 1),
     xlab = paste0("position across outlet channel (", um, ")"), ylab = "c", main = "outlet, x = 0")
lines(yo, step_prof(tree(c(0, 1, 0), 8), yo), type = "s", lwd = 3, col = oi["orange"])
legend("topleft", c("inlets (0, 1)", "inlets (0, 1, 0)"), col = c(oi["blue"], oi["orange"]),
       lwd = 3, bty = "n", cex = 1.2)
dev.off()
cat("fig11 outlet (0,1):", round(levels[[7]], 4), "\n")
cat("fig11 outlet (0,1,0):", round(tree(c(0, 1, 0), 8), 4), "\n")

## ---- fig12: outlet downstream --------------------------------------------------------------

Wo <- 400; Mo <- 160; dyo <- Wo / Mo; yo2 <- (1:Mo - 0.5) * dyo
uo <- pois(yo2, w = Wo, ubar = 200)
svg("fig12-outlet.svg", width = 13, height = 4.8)
par_big(mfrow = c(1, 2), mar = c(4.8, 5, 3, 1))
for (inl in list(c(0, 1), c(0, 1, 0))) {
  cs <- tree(inl, 8); co <- cs[ceiling(yo2 / 50)]
  plot(NA, xlim = c(0, Wo), ylim = c(0, 1), xlab = paste0("y (", um, ")"), ylab = "c",
       main = paste0("inlets (", paste(inl, collapse = ", "), ")"))
  xm <- c(0, 1000, 5000, 20000); cl <- c(oi["grey"], oi["sky"], oi["blue"], oi["red"])
  cc <- co; prev <- 0
  for (j in seq_along(xm)) {
    cc <- march(cc, uo, 100, dyo, 50, (xm[j] - prev) / 50); prev <- xm[j]
    lines(yo2, cc, lwd = 3, col = cl[j], type = if (j == 1) "s" else "l")
  }
  if (length(inl) == 2) legend("topleft", paste0("x = ", xm / 1000, " mm"), col = cl, lwd = 3, bty = "n", cex = 1.2)
}
dev.off()

## ---- fig13: fit D ---------------------------------------------------------------------------

profile <- function(D, M = 50, dx = 25, xobs = 1000, nbin = 25) {
  dy <- w / M; y <- (1:M - 0.5) * dy
  c <- march(ifelse(y < w / 2, 1, 0), pois(y), D, dy, dx, round(xobs / dx))
  colMeans(matrix(c, M / nbin))                     # average over 4 um pixels
}
cobs <- c(0.881, 0.906, 0.864, 0.859, 0.852, 0.828, 0.824, 0.758, 0.753, 0.681, 0.617,
          0.578, 0.496, 0.457, 0.383, 0.346, 0.269, 0.216, 0.153, 0.157, 0.112, 0.180,
          0.100, 0.089, 0.125)
yb <- (1:25 - 0.5) * 4
sse <- function(D) sum((cobs - profile(D))^2)
fit <- optimize(sse, c(10, 3000), tol = 1e-6)
cat("fig13 fit", fit$minimum, fit$objective, "\n")
Dg <- seq(100, 1200, by = 10); Sg <- sapply(Dg, sse)
svg("fig13-fit.svg", width = 13, height = 4.8)
par_big(mfrow = c(1, 2), mar = c(4.8, 5.5, 3, 1))
plot(yb, cobs, pch = 19, cex = 1.4, xlab = paste0("y (", um, ")"), ylab = "c (fluorescence, scaled)",
     main = "synthetic profile at x = 1 mm", ylim = c(0, 1))
lines(yb, profile(fit$minimum), lwd = 3, col = oi["blue"])
lines(yb, profile(150), lwd = 2, lty = 2, col = oi["orange"])
lines(yb, profile(1000), lwd = 2, lty = 3, col = oi["red"])
legend("topright", c(sprintf("fit, D = %.0f", fit$minimum), "D = 150", "D = 1000"),
       col = c(oi["blue"], oi["orange"], oi["red"]), lty = 1:3, lwd = c(3, 2, 2), bty = "n", cex = 1.2)
plot(Dg, Sg, type = "l", lwd = 3, log = "y", xlab = paste0("D (", um, "²/s)"), ylab = "SSE(D)",
     main = "one parameter: Brent's method")
abline(v = fit$minimum, lty = 2, col = oi["blue"])
dev.off()
