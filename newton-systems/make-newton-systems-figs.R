# Figures for the newton-systems/ unit ("Newton for systems of equations").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd newton-systems && Rscript make-newton-systems-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Running example: the symmetric two-gene toggle switch
#   du/dt = alpha / (1 + v^n) - u,   dv/dt = alpha / (1 + u^n) - v
# with alpha = 3, n = 2 (steady states (2.618, 0.382), (1.213, 1.213),
# (0.382, 2.618)). Makes:
#   fig01-toggle-circuit.svg   the circuit: two genes repressing each other
#   fig02-steady-states.svg    the curves F1 = 0 and F2 = 0 and their 3 crossings
#   fig03-partials.svg         partial derivatives as slopes of slices
#   fig04-newton-path.svg      Newton iterates from two starts
#   fig05-error-iter.svg       error vs iteration: exact and finite-difference J
#   fig06-basins.svg           which steady state Newton finds vs where the cell goes
#   fig07-fd-error.svg         error of the finite-difference Jacobian vs h
#   fig08-alpha-scan.svg       steady states vs promoter strength alpha

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

alpha <- 3; n <- 2
Fx <- function(x, a = alpha) c(a / (1 + x[2]^n) - x[1], a / (1 + x[1]^n) - x[2])
Jx <- function(x, a = alpha) {
  u <- x[1]; v <- x[2]
  matrix(c(-1, -a * n * v^(n - 1) / (1 + v^n)^2,
           -a * n * u^(n - 1) / (1 + u^n)^2, -1), 2, byrow = TRUE)
}
fdJ <- function(F, x, h = sqrt(.Machine$double.eps)) {
  f0 <- F(x); Jm <- matrix(0, length(f0), length(x))
  for (j in seq_along(x)) {
    xh <- x; hj <- h * max(abs(x[j]), 1); xh[j] <- x[j] + hj
    Jm[, j] <- (F(xh) - f0) / hj
  }
  Jm
}
newton_path <- function(x, Jf = Jx, maxit = 30, tol = 1e-14) {
  P <- matrix(x, nrow = 1)
  for (k in 1:maxit) {
    d <- solve(Jf(x), -Fx(x)); x <- x + d; P <- rbind(P, x)
    if (sqrt(sum(d^2)) < tol) break
  }
  P
}
x_hi  <- c((3 + sqrt(5)) / 2, (3 - sqrt(5)) / 2)
x_lo  <- rev(x_hi)
x_mid <- rep(uniroot(function(u) u^3 + u - 3, c(0, 3), tol = 1e-15)$root, 2)

big <- function() par(cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)

## ---- fig01: the circuit ---------------------------------------------------
svg("fig01-toggle-circuit.svg", width = 10, height = 3.6)
par(mar = c(0.2, 0.2, 0.2, 0.2))
plot(NA, xlim = c(0, 10), ylim = c(0, 3.6), axes = FALSE, xlab = "", ylab = "")
rect(0.6, 1.3, 3.4, 2.3, col = oi["sky"], border = NA)
rect(6.6, 1.3, 9.4, 2.3, col = oi["orange"], border = NA)
text(2, 1.8, "gene U", cex = 1.8); text(8, 1.8, "gene V", cex = 1.8)
# repression arrows: flat heads (angle = 90)
arrows(2.6, 2.5, 7.4, 2.5, angle = 90, length = 0.15, lwd = 4, col = oi["blue"])
arrows(7.4, 1.1, 2.6, 1.1, angle = 90, length = 0.15, lwd = 4, col = oi["red"])
text(5, 2.9, "protein U represses gene V", cex = 1.4, col = oi["blue"])
text(5, 0.65, "protein V represses gene U", cex = 1.4, col = oi["red"])
text(5, 0.15, "(Gardner, Cantor & Collins 2000: LacI and TetR in E. coli)", cex = 1.2, col = "grey30")
dev.off()

## ---- fig02: the two curves and the three steady states -------------------
nullcl <- function() {
  s <- seq(-0.5, 4.5, length.out = 400)
  lines(alpha / (1 + s^2), s, col = oi["blue"], lwd = 3)   # F1 = 0: u = 3/(1+v^2)
  lines(s, alpha / (1 + s^2), col = oi["red"], lwd = 3)    # F2 = 0: v = 3/(1+u^2)
}
svg("fig02-steady-states.svg", width = 7.5, height = 6.5)
par(mar = c(5, 5.2, 1, 1)); big()
plot(NA, xlim = c(0, 3.5), ylim = c(0, 3.5), xlab = "u (protein U)",
     ylab = "v (protein V)", asp = 1)
nullcl()
pts <- rbind(x_hi, x_mid, x_lo)
points(pts, pch = 21, bg = c(oi["green"], "white", oi["purple"]), cex = 2.6, lwd = 2)
text(x_hi[1], x_hi[2] + 0.3, "U high, V low", cex = 1.3)
text(x_lo[1] + 0.2, x_lo[2] + 0.05, "V high, U low", cex = 1.3, adj = 0)
text(x_mid[1] + 0.2, x_mid[2] + 0.2, "both medium", cex = 1.3, adj = 0)
legend("topright", c(expression(F[1] == 0 ~ ": " ~ u == 3/(1+v^2)),
                     expression(F[2] == 0 ~ ": " ~ v == 3/(1+u^2))),
       col = c(oi["blue"], oi["red"]), lwd = 3, cex = 1.2, bty = "n")
dev.off()

## ---- fig03: partial derivatives as slopes of slices -----------------------
svg("fig03-partials.svg", width = 12, height = 4.8)
par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1)); big()
u <- seq(1, 5, length.out = 200)
plot(u, alpha / (1 + 0.5^2) - u, type = "l", lwd = 3, col = oi["blue"],
     xlab = "u", ylab = expression(F[1](u, 0.5)),
     main = "hold v = 0.5, move u")
abline(a = alpha / 1.25, b = -1, lty = 2, lwd = 2, col = oi["orange"])
points(3, -0.6, pch = 19, cex = 1.8)
text(3.3, 0.3, expression(partialdiff * F[1] / partialdiff * u == -1), cex = 1.4, adj = 0)
v <- seq(0, 2, length.out = 200)
plot(v, alpha / (1 + v^2) - 3, type = "l", lwd = 3, col = oi["blue"],
     xlab = "v", ylab = expression(F[1](3, v)),
     main = "hold u = 3, move v")
abline(a = -0.6 + 1.92 * 0.5, b = -1.92, lty = 2, lwd = 2, col = oi["orange"])
points(0.5, -0.6, pch = 19, cex = 1.8)
text(0.75, -0.25, expression(partialdiff * F[1] / partialdiff * v == -1.92), cex = 1.4, adj = 0)
dev.off()

## ---- fig04: Newton paths --------------------------------------------------
svg("fig04-newton-path.svg", width = 8, height = 6.5)
par(mar = c(5, 5.2, 1, 1)); big()
plot(NA, xlim = c(0, 5.2), ylim = c(-1, 3.2), xlab = "u (protein U)",
     ylab = "v (protein V)", asp = 1)
abline(h = 0, v = 0, col = "grey80")
nullcl()
points(pts, pch = 21, bg = c(oi["green"], "white", oi["purple"]), cex = 2.2, lwd = 2)
P1 <- newton_path(c(3, 0.5)); P2 <- newton_path(c(2, 0.5))
for (P in list(list(P1, oi["orange"]), list(P2, "black"))) {
  M <- P[[1]]; k <- nrow(M)
  long <- sqrt(rowSums((M[-1, , drop = FALSE] - M[-k, , drop = FALSE])^2)) > 0.03
  arrows(M[-k, 1][long], M[-k, 2][long], M[-1, 1][long], M[-1, 2][long],
         length = 0.1, lwd = 2, col = P[[2]])
  points(M[1, 1], M[1, 2], pch = 15, cex = 1.8, col = P[[2]])
}
text(3.05, 0.75, "start (3, 0.5)", cex = 1.3, col = oi["orange"], adj = 0)
text(1.9, 0.3, "start (2, 0.5)", cex = 1.3, adj = 1)
text(3.6, -0.9, "v < 0: not a concentration", cex = 1.2, col = "grey30")
dev.off()

## ---- fig05: error vs iteration --------------------------------------------
svg("fig05-error-iter.svg", width = 8, height = 5.5)
par(mar = c(5, 5.5, 1, 1)); big()
e_ex <- apply(newton_path(c(3, 0.5)), 1, function(x) sqrt(sum((x - x_hi)^2)))
e_fd <- apply(newton_path(c(3, 0.5), Jf = function(x) fdJ(Fx, x)), 1,
              function(x) sqrt(sum((x - x_hi)^2)))
e_ex[e_ex == 0] <- 1e-17; e_fd[e_fd == 0] <- 1e-17
plot(seq_along(e_ex) - 1, e_ex, log = "y", type = "b", pch = 19, lwd = 2,
     col = oi["blue"], ylim = c(1e-17, 1), xlim = c(0, 6), yaxt = "n",
     xlab = "iteration k", ylab = expression("error " * "||" * x[k] - x^"*" * "||"))
axis(2, at = 10^c(-16, -12, -8, -4, 0), labels = expression(10^-16, 10^-12, 10^-8, 10^-4, 1), las = 1)
lines(seq_along(e_fd) - 1, e_fd, type = "b", pch = 17, lwd = 2, lty = 2, col = oi["orange"])
abline(h = .Machine$double.eps, col = oi["grey"], lty = 3)
legend("topright", c("exact Jacobian", "finite-difference Jacobian"),
       col = c(oi["blue"], oi["orange"]), pch = c(19, 17), lty = c(1, 2), lwd = 2,
       cex = 1.2, bty = "n")
dev.off()

## ---- fig06: basins: Newton vs the ODE --------------------------------------
g <- seq(0.005, 4, length.out = 300)
U <- matrix(rep(g, times = length(g)), length(g))         # U[i, j] = g[i]
V <- matrix(rep(g, each = length(g)), length(g))          # V[i, j] = g[j]
u <- U; v <- V
for (k in 1:60) {                     # vectorised Newton with a 2x2 Cramer solve
  f1 <- alpha / (1 + v^2) - u; f2 <- alpha / (1 + u^2) - v
  b <- -2 * alpha * v / (1 + v^2)^2; c <- -2 * alpha * u / (1 + u^2)^2
  det <- 1 - b * c                    # J = [-1 b; c -1]
  du <- (f1 + b * f2) / det           # solve J d = -F
  dv <- (c * f1 + f2) / det
  u <- u + du; v <- v + dv
}
cls <- function(u, v) {
  out <- matrix(4, nrow(u), ncol(u))  # 4 = no convergence
  ok <- is.finite(u) & is.finite(v)
  out[ok & abs(u - x_hi[1]) + abs(v - x_hi[2]) < 1e-6] <- 1
  out[ok & abs(u - x_mid[1]) + abs(v - x_mid[2]) < 1e-6] <- 2
  out[ok & abs(u - x_lo[1]) + abs(v - x_lo[2]) < 1e-6] <- 3
  out
}
Bn <- cls(u, v)
u <- U; v <- V
for (k in 1:2000) {                   # Euler steps of the ODE, dt = 0.05, to t = 100
  du <- alpha / (1 + v^2) - u; dv <- alpha / (1 + u^2) - v
  u <- u + 0.05 * du; v <- v + 0.05 * dv
}
Bo <- cls(u, v)
cols <- c(oi["green"], "grey85", oi["purple"], "black")
svg("fig06-basins.svg", width = 12, height = 5.8)
par(mfrow = c(1, 2), mar = c(5, 5.2, 3, 1)); big()
for (B in list(list(Bn, "Newton converges to ..."), list(Bo, "the cell (ODE) ends in ..."))) {
  image(g, g, B[[1]], col = cols, breaks = 0.5:4.5, xlab = "start u", ylab = "start v",
        main = B[[2]], asp = 1, useRaster = TRUE)
  points(pts, pch = 21, bg = "white", cex = 2, lwd = 2)
}
legend("topright", c("U high", "middle", "V high"), fill = cols[1:3], cex = 1.2, bg = "white")
dev.off()
cat("Newton basin shares (hi, mid, lo, fail):", round(tabulate(Bn, 4) / length(Bn), 3), "\n")

## ---- fig07: finite-difference Jacobian error vs h --------------------------
hs <- 10^seq(-1, -14, by = -0.25)
err <- sapply(hs, function(h) max(abs(fdJ(Fx, c(3, 0.5), h) - Jx(c(3, 0.5)))))
svg("fig07-fd-error.svg", width = 8, height = 5.5)
par(mar = c(5, 5.5, 1, 1)); big()
plot(hs, err, log = "xy", type = "b", pch = 19, lwd = 2, col = oi["blue"],
     xlab = "step h", ylab = "max error in J", xaxt = "n", yaxt = "n")
axis(1, at = 10^c(-14, -11, -8, -5, -2), labels = expression(10^-14, 10^-11, 10^-8, 10^-5, 10^-2))
axis(2, at = 10^c(-8, -6, -4, -2), labels = expression(10^-8, 10^-6, 10^-4, 10^-2), las = 1)
abline(v = sqrt(.Machine$double.eps), lty = 2, col = oi["orange"], lwd = 2)
text(sqrt(.Machine$double.eps) * 1.6, 1e-2, expression(h == sqrt(epsilon)), adj = 0,
     cex = 1.3, col = oi["orange"])
text(3e-3, 3e-6, "truncation\n(slope 1)", cex = 1.2)
text(1e-12, 3e-6, "rounding\n(slope -1)", cex = 1.2)
dev.off()

## ---- fig08: steady states vs alpha -----------------------------------------
svg("fig08-alpha-scan.svg", width = 8, height = 5.5)
par(mar = c(5, 5.5, 1, 1)); big()
plot(NA, xlim = c(1, 5), ylim = c(0, 5), xlab = expression("promoter strength " * alpha),
     ylab = "u at steady state")
starts <- as.matrix(expand.grid(seq(0.1, 5, length.out = 7), seq(0.1, 5, length.out = 7)))
for (a in seq(1, 5, by = 0.05)) {
  found <- NULL
  for (i in seq_len(nrow(starts))) {
    x <- starts[i, ]
    for (k in 1:50) {
      d <- tryCatch(solve(Jx(x, a), -Fx(x, a)), error = function(e) c(NaN, NaN))
      x <- x + d
      if (!all(is.finite(x))) break
      if (sqrt(sum(d^2)) < 1e-12) break
    }
    if (all(is.finite(x)) && sqrt(sum(Fx(x, a)^2)) < 1e-10) found <- rbind(found, x)
  }
  found <- unique(round(found, 6))
  sym <- abs(found[, 1] - found[, 2]) < 1e-6
  points(rep(a, sum(sym)), found[sym, 1], pch = 19, cex = 0.9,
         col = ifelse(a > 2, "grey50", "black"))
  points(rep(a, sum(!sym)), found[!sym, 1], pch = 19, cex = 0.9,
         col = ifelse(found[!sym, 1] > found[!sym, 2], oi["green"], oi["purple"]))
}
abline(v = 2, lty = 2)
text(2.05, 4.7, expression(alpha == 2), adj = 0, cex = 1.3)
text(1.5, 4.0, "one steady state", cex = 1.3)
text(2.9, 4.6, "three: bistable", cex = 1.3)
dev.off()
