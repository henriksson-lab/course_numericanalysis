# Figures for the reaction-diffusion/ unit ("Reaction-diffusion and pattern formation").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd reaction-diffusion && Rscript make-reaction-diffusion-figs.R
#
# Base R graphics only (grDevices::svg) plus Matrix::expm for the exact
# sub-steps in fig04. Colours are Okabe-Ito (colour-blind safe). Takes about
# half a minute (the 2D runs). Makes:
#   fig01-bicoid.svg       embryo with an anterior source, and the exponential
#                          gradient with its decay length lambda
#   fig02-steady.svg       steady gradient on the grid vs the exact cosh formula,
#                          and error vs dx: ghost point (order 2) vs one-sided
#                          boundary (order 1)
#   fig03-formation.svg    the gradient forming (IMEX, dt = 1 min) and the
#                          time course at three positions
#   fig04-splitting.svg    error vs dt at t = 60 min: Lie, Strang, IMEX
#   fig05-tstep.svg        steady gradient with the temperature step (assumed
#                          rate changes) and where it crosses a threshold
#   fig06-dispersion.svg   largest growth rate vs wavenumber q (Schnakenberg)
#                          for several diffusion ratios d
#   fig07-turing.svg       the 2D pattern forming (u at t = 10, 30, 200)
#   fig08-artefacts.svg    converged pattern, coarse grid, and a "pattern"
#                          made by a too-large time step (d = 5: no Turing)
#
# All numbers are computed here with the same functions as on the slides.

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- 1D synthesis-diffusion-degradation (units: µm, min) ----------------------

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

D <- 100; k <- 0.01; L <- 500          # µm^2/min, 1/min, µm
lam <- sqrt(D / k); j <- D / lam        # lambda = 100 µm; c(0) close to 1
exact <- function(x) j * lam / D * cosh((L - x) / lam) / sinh(L / lam)

# nodes x_0..x_N, ghost points at both ends; returns the positive operator
# -D d2/dx2 as three diagonals and the boundary source vector
grid1d <- function(N) {
  dx <- L / N; r <- D / dx^2
  list(N = N, dx = dx, x = (0:N) * dx,
       lo = c(0, rep(-r, N - 1), -2 * r), di = rep(2 * r, N + 1),
       up = c(-2 * r, rep(-r, N - 1), 0), b = c(2 * j / dx, rep(0, N)))
}
steady <- function(N, kx = k, onesided = FALSE) {
  g <- grid1d(N); kx <- rep(kx, length.out = N + 1)
  di <- g$di + kx; up <- g$up; rhs <- g$b
  if (onesided) {                    # (c1 - c0)/dx = -j/D, first order
    di[1] <- 1; up[1] <- -1; rhs[1] <- j * g$dx / D
  }
  list(x = g$x, c = thomas(g$lo, di, up, rhs))
}
imex <- function(c, g, dt, nsteps, f) {
  for (n in seq_len(nsteps))
    c <- thomas(dt * g$lo, 1 + dt * g$di, dt * g$up, c + dt * (f(c) + g$b))
  c
}

## ---- fig01: Bicoid set-up --------------------------------------------------------

svg("fig01-bicoid.svg", width = 12, height = 5)
layout(matrix(1:2, 1), widths = c(1, 1.15))
par(mar = c(1, 1, 1, 1))
plot(NA, xlim = c(-30, 530), ylim = c(-150, 150), asp = 1, axes = FALSE, xlab = "", ylab = "")
th <- seq(0, 2 * pi, length.out = 300)
xs <- seq(0, 500, by = 2)
cols <- colorRampPalette(c("white", oi["blue"]))(101)
for (xx in xs) {
  hgt <- 95 * sqrt(pmax(0, 1 - ((xx - 250) / 250)^2))
  segments(xx, -hgt, xx, hgt, col = cols[1 + round(100 * exact(xx))], lwd = 2.2)
}
lines(250 + 250 * cos(th), 95 * sin(th), lwd = 3)
points(8, 0, pch = 19, cex = 2.5, col = oi["red"])
text(15, 120, "anterior: bcd mRNA (source)", pos = 4, cex = 1.3, col = oi["red"])
text(500, -125, "posterior", pos = 2, cex = 1.3)
text(250, -142, "embryo ≈ 500 µm", cex = 1.3)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
xx <- seq(0, 500, by = 2)
plot(xx, exact(xx), type = "l", lwd = 4, col = oi["blue"], las = 1,
     xlab = "distance from the anterior pole x (µm)", ylab = "Bicoid c / c(0)")
abline(h = exp(-1), lty = 3); segments(lam, 0, lam, exp(-1), lty = 3)
text(lam + 5, 0.45, expression(lambda == sqrt(D / k) ~ "=" ~ 100 ~ mu * m), pos = 4, cex = 1.4)
text(250, 0.85, "steady state of\n dc/dt = D c'' - k c", cex = 1.3, pos = 4)
dev.off()

## ---- fig02: steady state on the grid, and its order -------------------------------

s50 <- steady(50); s5 <- steady(5)
cat("steady N=50 max err", max(abs(s50$c - exact(s50$x))), "\n")
Ns <- c(10, 20, 40, 80, 160, 320, 640)
eg <- sapply(Ns, function(N) { s <- steady(N); max(abs(s$c - exact(s$x))) })
eo <- sapply(Ns, function(N) { s <- steady(N, onesided = TRUE); max(abs(s$c - exact(s$x))) })
print(rbind(dx = L / Ns, ghost = signif(eg, 3), onesided = signif(eo, 3)))
svg("fig02-steady.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(xx, exact(xx), type = "l", lwd = 4, col = "grey60", las = 1, ylim = c(0, 1.05),
     xlab = "x (µm)", ylab = "steady c")
points(s50$x, s50$c, pch = 19, col = oi["blue"], cex = 1.1)
points(s5$x, s5$c, pch = 17, col = oi["red"], cex = 1.8)
legend("topright", c("exact (cosh formula)", "grid, dx = 10 µm", "grid, dx = 100 µm = lambda"),
       lwd = c(4, NA, NA), pch = c(NA, 19, 17), col = c("grey60", oi["blue"], oi["red"]),
       cex = 1.2, bty = "n")
plot(L / Ns, eg, log = "xy", type = "b", pch = 19, lwd = 3, col = oi["blue"], las = 1,
     ylim = range(c(eg, eo)), xlab = "grid spacing dx (µm)", ylab = "max error")
lines(L / Ns, eo, type = "b", pch = 17, lwd = 3, col = oi["red"])
text(L / Ns[2], eg[2], "ghost point: slope 2", pos = 4, cex = 1.25, col = oi["blue"])
text(L / Ns[5], eo[5], "one-sided boundary: slope 1", pos = 3, cex = 1.25, col = oi["red"])
dev.off()

## ---- fig03: formation -------------------------------------------------------------

g <- grid1d(50)
fk <- function(c) -k * c
tt <- 0:600; ix <- c(1, 11, 21)
cs <- steady(50)$c
traj <- matrix(NA, length(tt), 3); c <- rep(0, 51); traj[1, ] <- 0
prof <- list()
for (s in 2:length(tt)) {
  c <- imex(c, g, 1, 1, fk); traj[s, ] <- c[ix] / cs[ix]
  if (tt[s] %in% c(10, 30, 100, 300)) prof[[as.character(tt[s])]] <- c
}
tau_num <- colSums(1 - traj) * 1                  # rectangle rule, dt = 1 min
cat("accumulation times (num, BE dt=1):", round(tau_num, 1), " formula:", (1 + c(0, 100, 200) / lam) / (2 * k), "\n")
svg("fig03-formation.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(g$x, cs, type = "l", lwd = 4, col = "black", las = 1, ylim = c(0, 1.05),
     xlab = "x (µm)", ylab = "c")
pc <- oi[c("sky", "green", "orange", "red")]
for (i in 1:4) lines(g$x, prof[[i]], lwd = 3, col = pc[i])
legend("topright", c(paste("t =", names(prof), "min"), "steady state"), lwd = c(3, 3, 3, 3, 4),
       col = c(pc, "black"), cex = 1.2, bty = "n")
plot(tt, traj[, 1], type = "l", lwd = 4, col = oi["blue"], ylim = c(0, 1), las = 1,
     xlab = "time (min)", ylab = "c / steady c")
lines(tt, traj[, 2], lwd = 4, col = oi["orange"]); lines(tt, traj[, 3], lwd = 4, col = oi["red"])
legend("bottomright", c("x = 0", "x = 100 µm", "x = 200 µm"), lwd = 4,
       col = oi[c("blue", "orange", "red")], cex = 1.2, bty = "n")
dev.off()

## ---- fig04: splitting error -----------------------------------------------------------
## A = diffusion matrix with ghost points (incl. the boundary flux, as an affine
## flow done exactly with expm of an augmented matrix); reaction -k(x) c done
## exactly. Reference: expm of the full linear system.

library(Matrix)
x <- g$x; m <- 51
A <- matrix(0, m, m); A[cbind(1:m, 1:m)] <- -g$di
A[cbind(2:m, 1:(m - 1))] <- -g$lo[-1]; A[cbind(1:(m - 1), 2:m)] <- -g$up[-m]
kx <- ifelse(x < L / 2, 0.01, 0.02)
M <- A - diag(kx); csx <- -solve(M, g$b)
tend <- 60
cref <- csx + as.vector(as.matrix(expm(Matrix(M * tend))) %*% (0 - csx))
Aug <- rbind(cbind(A, g$b), 0)
hs <- c(20, 10, 5, 2.5, 1.25, 0.625)
errs <- t(sapply(hs, function(h) {
  E <- as.matrix(expm(Matrix(Aug * h))); dif <- function(c) as.vector(E %*% c(c, 1))[1:m]
  c1 <- c2 <- c3 <- rep(0, m)
  for (s in seq_len(tend / h)) {
    c1 <- exp(-kx * h) * dif(c1)                         # Lie
    c2 <- exp(-kx * h / 2) * dif(exp(-kx * h / 2) * c2)  # Strang
    c3 <- imex(c3, g, h, 1, function(c) -kx * c)         # IMEX
  }
  c(max(abs(c1 - cref)), max(abs(c2 - cref)), max(abs(c3 - cref)))
}))
print(cbind(dt = hs, signif(errs, 3)))
svg("fig04-splitting.svg", width = 8, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
matplot(hs, errs, log = "xy", type = "b", pch = c(15, 19, 17), lty = 1, lwd = 3,
        col = oi[c("red", "blue", "green")], las = 1,
        xlab = "time step dt (min)", ylab = "max error at t = 60 min")
legend("bottomright", c("Lie splitting (slope 1)", "Strang splitting (slope 2)", "IMEX Euler (slope 1)"),
       pch = c(15, 19, 17), lwd = 3, col = oi[c("red", "blue", "green")], cex = 1.2, bty = "n")
dev.off()

## ---- fig05: temperature step ----------------------------------------------------------
## Assumed: in the warm (27 °C) half, k is 2x and D is 1.2x their 20 °C values.

Dloc <- function(N, Dx, kx) {           # variable D: flux form, ghost points
  dx <- L / N; xn <- (0:N) * dx; Dm <- Dx((xn[-1] + xn[-(N + 1)]) / 2)   # D at midpoints
  lo <- c(0, -Dm[-N] / dx^2, -2 * Dm[N] / dx^2)
  up <- c(-2 * Dm[1] / dx^2, -Dm[-1] / dx^2, 0)
  di <- c(2 * Dm[1], Dm[-N] + Dm[-1], 2 * Dm[N]) / dx^2 + kx(xn)
  rhs <- c(2 * j / dx, rep(0, N))
  list(x = xn, c = thomas(lo, di, up, rhs))
}
N5 <- 500
cold <- function(x) rep(1, length(x)); warm <- function(x) rep(1, length(x))
s_cold <- Dloc(N5, function(x) D + 0 * x, function(x) k + 0 * x)
s_warm <- Dloc(N5, function(x) 1.2 * D + 0 * x, function(x) 2 * k + 0 * x)
s_pw <- Dloc(N5, function(x) ifelse(x < L / 2, D, 1.2 * D), function(x) ifelse(x < L / 2, k, 2 * k))
s_aw <- Dloc(N5, function(x) ifelse(x < L / 2, 1.2 * D, D), function(x) ifelse(x < L / 2, 2 * k, k))
xT <- function(s, thr = 0.1) approx(s$c, s$x, thr)$y
cat("threshold 0.1 at: cold", xT(s_cold), "warm", xT(s_warm), "post warm", xT(s_pw), "ant warm", xT(s_aw), "\n")
svg("fig05-tstep.svg", width = 10, height = 5.2)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(s_cold$x, s_cold$c, type = "l", lwd = 4, col = oi["blue"], log = "y", ylim = c(0.005, 1.1),
     las = 1, xlab = "x (µm)", ylab = "steady c (log scale)")
lines(s_warm$x, s_warm$c, lwd = 4, col = oi["red"])
lines(s_pw$x, s_pw$c, lwd = 3, lty = 2, col = oi["purple"])
lines(s_aw$x, s_aw$c, lwd = 3, lty = 2, col = oi["orange"])
abline(h = 0.1, lty = 3); abline(v = L / 2, col = "grey70")
legend("bottomleft", c("all 20 °C", "all 27 °C", "posterior half warm", "anterior half warm"),
       lwd = c(4, 4, 3, 3), lty = c(1, 1, 2, 2), col = oi[c("blue", "red", "purple", "orange")],
       cex = 1.2, bty = "n")
text(300, 0.115, "threshold c = 0.1", pos = 3, cex = 1.2)
dev.off()

## ---- Schnakenberg: dispersion relation -------------------------------------------

a <- 0.1; b <- 0.9; us <- a + b; vs <- b / us^2
J <- matrix(c(-1 + 2 * us * vs, -2 * us * vs, us^2, -us^2), 2)
growth <- function(q, d) sapply(q, function(qq) max(Re(eigen(J - qq^2 * diag(c(1, d)))$values)))
q <- seq(0, 1.2, by = 0.002)
svg("fig06-dispersion.svg", width = 9, height = 5.2)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
ds <- c(1, 5, 10, 20); dc <- oi[c("grey", "sky", "orange", "red")]
plot(q, growth(q, 1), type = "l", lwd = 4, col = dc[1], ylim = c(-0.6, 0.3), las = 1,
     xlab = "wavenumber q (2 pi / wavelength)", ylab = expression("largest Re " * lambda(q)))
for (i in 2:4) lines(q, growth(q, ds[i]), lwd = 4, col = dc[i])
abline(h = 0, lty = 2)
g20 <- growth(q, 20); qm <- q[which.max(g20)]
cat("d = 20: fastest q", qm, "wavelength", 2 * pi / qm, "rate", max(g20), "band", range(q[g20 > 0]), "\n")
points(qm, max(g20), pch = 19, cex = 1.5, col = oi["red"])
legend("bottomleft", paste("d =", ds), lwd = 4, col = dc, cex = 1.2, bty = "n")
dev.off()

## ---- 2D Turing simulation (IMEX, FFT for the implicit 5-point diffusion) -----------

lcg <- function(n, seed = 1) {
  z <- numeric(n); s <- seed
  for (i in seq_len(n)) { s <- (16807 * s) %% 2147483647; z[i] <- s / 2147483647 }
  z
}
schnak <- function(N, dt, tend, d = 20, Lbox = 80) {
  dx <- Lbox / N; s <- 4 / dx^2 * sin(pi * (0:(N - 1)) / N)^2
  kap <- outer(s, s, "+")
  z <- lcg(2 * N^2) - 0.5
  u <- matrix(1 + 0.02 * z[1:N^2], N); v <- matrix(0.9 + 0.02 * z[-(1:N^2)], N)
  sd <- function(w, Dw) Re(fft(fft(w) / (1 + dt * Dw * kap), inverse = TRUE)) / N^2
  for (n in seq_len(round(tend / dt))) {
    f <- a - u + u^2 * v; g <- b - u^2 * v
    u <- sd(u + dt * f, 1); v <- sd(v + dt * g, d)
  }
  u
}
spots <- function(u) {
  N <- nrow(u); ismax <- u > mean(u)
  for (i in -1:1) for (jj in -1:1) if (i != 0 || jj != 0)
    ismax <- ismax & u > u[(0:(N - 1) + i) %% N + 1, (0:(N - 1) + jj) %% N + 1]
  sum(ismax)
}
pal <- hcl.colors(64, "viridis")
img <- function(u, main, zl = c(0, 3.3)) {
  N <- nrow(u); xx <- (1:N - 0.5) * 80 / N
  image(xx, xx, pmin(pmax(u, zl[1]), zl[2]), col = pal, zlim = zl, asp = 1, axes = FALSE,
        xlab = "", ylab = "", main = main)
  box()
}
svg("fig07-turing.svg", width = 12, height = 4.4)
par(mfrow = c(1, 3), mar = c(1, 1, 3.5, 1), cex.main = 1.8)
for (te in c(10, 30, 200)) {
  u <- schnak(80, 0.1, te)
  img(u, sprintf("t = %d  (u from %.2f to %.2f)", te, min(u), max(u)), zl = range(u))
}
dev.off()

u_ref <- schnak(80, 0.1, 200); u_co <- schnak(40, 0.1, 200); u_fi <- schnak(160, 0.1, 200)
u_d5 <- schnak(80, 0.4, 150, d = 5); u_d5ok <- schnak(80, 0.1, 150, d = 5)
cat("spots: dx=1", spots(u_ref), " dx=2", spots(u_co), " dx=0.5", spots(u_fi), "\n")
cat("d=5 dt=0.4 range", range(u_d5), "  d=5 dt=0.1 range", range(u_d5ok), "\n")
svg("fig08-artefacts.svg", width = 12, height = 4.4)
par(mfrow = c(1, 3), mar = c(1, 1, 3.5, 1), cex.main = 1.7)
img(u_ref, sprintf("dx = 1, dt = 0.1: %d spots", spots(u_ref)))
img(u_co, sprintf("dx = 2: %d spots", spots(u_co)))
img(u_d5, "d = 5 (no Turing), dt = 0.4", zl = range(u_d5))
dev.off()
