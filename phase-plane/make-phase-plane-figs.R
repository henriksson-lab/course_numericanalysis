# Figures for the phase-plane/ unit ("Equilibria, stability and the phase plane").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd phase-plane && Rscript make-phase-plane-figs.R
#
# Base R graphics only (grDevices::svg) plus deSolve for trajectories.
# Colours are Okabe-Ito (colour-blind safe).
# Running examples:
#   toggle switch  du/dt = alpha/(1+v^2) - u, dv/dt = alpha/(1+u^2) - v, alpha = 3
#                  (defined in newton-systems; steady states (2.618, 0.382),
#                  (1.213, 1.213), (0.382, 2.618))
#   biased toggle  as above with alpha_u = 3, alpha_v = 2.9
#   Sel'kov glycolysis  du/dt = -u + a v + u^2 v, dv/dt = b - a v - u^2 v,
#                  a = 0.08, b = 0.6 (limit cycle)
# Makes:
#   fig01-vector-field.svg   toggle: arrows of the flow, trajectories, steady states
#   fig02-nullclines.svg     toggle: nullclines, flow direction on them and between them
#   fig03-linear-zoo.svg     node, saddle, spiral, centre for linear 2x2 systems
#   fig04-trace-det.svg      the trace-determinant plane with our steady states
#   fig05-saddle.svg         zoom on the middle state: eigenvectors and trajectories
#   fig06-separatrix.svg     biased toggle: Newton basins vs ODE basins + separatrix
#   fig07-saddle-slow.svg    time courses from starts on / near the stable manifold
#   fig08-selkov.svg         Sel'kov: phase plane with limit cycle, and u(t)
#   fig09-bifurcation.svg    toggle: steady states and eigenvalues vs alpha

library(deSolve)

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")
big <- function() par(cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)

alpha <- 3
tog <- function(t, x, p, au = alpha, av = alpha)
  list(c(au / (1 + x[2]^2) - x[1], av / (1 + x[1]^2) - x[2]))
Jt <- function(x, au = alpha, av = alpha) {
  u <- x[1]; v <- x[2]
  matrix(c(-1, -2 * au * v / (1 + v^2)^2, -2 * av * u / (1 + u^2)^2, -1), 2, byrow = TRUE)
}
x_hi  <- c((3 + sqrt(5)) / 2, (3 - sqrt(5)) / 2)
x_lo  <- rev(x_hi)
x_mid <- rep(uniroot(function(u) u^3 + u - 3, c(0, 3), tol = 1e-15)$root, 2)

# normalised arrows of a vector field f(u, v) -> list(du, dv) on a grid
field <- function(f, ur, vr, nx = 15, len = 0.11, col = "grey55") {
  g <- expand.grid(u = seq(ur[1], ur[2], length.out = nx),
                   v = seq(vr[1], vr[2], length.out = nx))
  d <- f(g$u, g$v); s <- sqrt(d[[1]]^2 + d[[2]]^2); s[s == 0] <- 1
  sc <- len * diff(ur)
  arrows(g$u, g$v, g$u + sc * d[[1]] / s, g$v + sc * d[[2]] / s,
         length = 0.05, col = col, lwd = 1.2)
}
ftog <- function(u, v) list(alpha / (1 + v^2) - u, alpha / (1 + u^2) - v)
nullcl <- function(au = alpha, av = alpha, lwd = 3) {
  s <- seq(0, 5, length.out = 400)
  lines(au / (1 + s^2), s, col = oi["blue"], lwd = lwd)   # du/dt = 0
  lines(s, av / (1 + s^2), col = oi["red"], lwd = lwd)    # dv/dt = 0
}
states <- function(cex = 2.4) {
  points(rbind(x_hi, x_lo), pch = 21, bg = "black", cex = cex, lwd = 2)
  points(x_mid[1], x_mid[2], pch = 21, bg = "white", cex = cex, lwd = 2)
}
traj <- function(func, y0, tmax, col, lwd = 2.5, dt = 0.02, ...) {
  o <- ode(y0, seq(0, tmax, by = dt), func, NULL, ...)
  lines(o[, 2], o[, 3], col = col, lwd = lwd)
  k <- nrow(o) %/% 6
  arrows(o[k, 2], o[k, 3], o[k + 1, 2], o[k + 1, 3], length = 0.12, col = col, lwd = lwd)
  points(y0[1], y0[2], pch = 15, col = col, cex = 1.4)
  invisible(o)
}

## ---- fig01: vector field and trajectories ---------------------------------
svg("fig01-vector-field.svg", width = 7.5, height = 6.5)
par(mar = c(5, 5.2, 1, 1)); big()
plot(NA, xlim = c(0, 3.5), ylim = c(0, 3.5), xlab = "u (protein U)",
     ylab = "v (protein V)", asp = 1)
field(ftog, c(0.1, 3.4), c(0.1, 3.4))
starts <- list(c(0.5, 0.3), c(0.3, 0.5), c(3.3, 3.0), c(3.0, 3.3), c(0.1, 1.8),
               c(1.8, 0.1), c(2.0, 2.0))
cl <- c(oi["orange"], oi["sky"], oi["orange"], oi["sky"], oi["sky"], oi["orange"], oi["grey"])
for (i in seq_along(starts)) traj(tog, starts[[i]], 30, cl[i])
states()
dev.off()

## ---- fig02: nullclines ----------------------------------------------------
svg("fig02-nullclines.svg", width = 7.5, height = 6.5)
par(mar = c(5, 5.2, 1, 1)); big()
plot(NA, xlim = c(0, 3.5), ylim = c(0, 3.5), xlab = "u (protein U)",
     ylab = "v (protein V)", asp = 1)
nullcl()
# on the u-nullcline the flow is vertical, on the v-nullcline horizontal
for (vv in c(0.15, 0.8, 1.7, 2.4, 3.2)) {
  uu <- alpha / (1 + vv^2); d <- sign(alpha / (1 + uu^2) - vv)
  arrows(uu, vv, uu, vv + 0.22 * d, length = 0.09, lwd = 2.5, col = oi["blue"])
}
for (uu in c(0.15, 0.8, 1.7, 2.4, 3.2)) {
  vv <- alpha / (1 + uu^2); d <- sign(alpha / (1 + vv^2) - uu)
  arrows(uu, vv, uu + 0.22 * d, vv, length = 0.09, lwd = 2.5, col = oi["red"])
}
states()
text(0.45, 0.45, "u ↑  v ↑", cex = 1.3)
text(2.7, 2.9, "u ↓  v ↓", cex = 1.3)
# the two thin lenses between the curves
arrows(2.45, 1.5, 1.95, 0.8, length = 0.08, lwd = 1.5)
text(2.5, 1.62, "u \u2191  v \u2193", cex = 1.3)
arrows(1.45, 2.45, 0.75, 1.95, length = 0.08, lwd = 1.5)
text(1.55, 2.55, "u \u2193  v \u2191", cex = 1.3)
legend("topright", c("du/dt = 0  (u-nullcline)", "dv/dt = 0  (v-nullcline)"),
       col = c(oi["blue"], oi["red"]), lwd = 3, cex = 1.15, bg = "white")
dev.off()

## ---- fig03: the zoo of linear systems -------------------------------------
lin <- function(A) function(t, x, p) list(as.vector(A %*% x))
zoo <- list(
  list(A = matrix(c(-2, 0, 0, -0.5), 2), main = "stable node",
       sub = expression(lambda == -2 * "," ~ -0.5)),
  list(A = matrix(c(1, 0, 0, -1), 2), main = "saddle",
       sub = expression(lambda == 1 * "," ~ -1)),
  list(A = matrix(c(-0.3, -1, 1, -0.3), 2), main = "stable spiral",
       sub = expression(lambda == -0.3 %+-% i)),
  list(A = matrix(c(0, -1, 1, 0), 2), main = "centre",
       sub = expression(lambda == 0 %+-% i)))
svg("fig03-linear-zoo.svg", width = 13, height = 3.9)
par(mfrow = c(1, 4), mar = c(3.2, 1, 3, 1)); big()
for (z in zoo) {
  plot(NA, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, axes = FALSE, xlab = "", ylab = "",
       main = z$main)
  box(col = "grey70"); title(sub = z$sub, line = 1.8, cex.sub = 1.6)
  field(function(u, v) list(z$A[1, 1] * u + z$A[1, 2] * v, z$A[2, 1] * u + z$A[2, 2] * v),
        c(-0.95, 0.95), c(-0.95, 0.95), nx = 9, len = 0.06, col = "grey70")
  ang <- seq(0, 2 * pi, length.out = 13)[-13] + 0.2
  for (k in seq_along(ang)) {
    y0 <- 0.95 * c(cos(ang[k]), sin(ang[k]))
    if (z$main == "saddle") y0 <- c(c(-0.9, -0.3, 0.3, 0.9)[(k - 1) %% 4 + 1], ifelse(k > 6, 0.95, -0.95))
    if (z$main == "centre") y0 <- c(0.1 + 0.2 * k, 0)
    if (z$main %in% c("stable spiral", "centre") && k > 4) next
    o <- ode(y0, seq(0, ifelse(z$main == "centre", 6.3, 12), by = 0.02), lin(z$A), NULL)
    keep <- abs(o[, 2]) <= 1 & abs(o[, 3]) <= 1
    lines(o[keep, 2], o[keep, 3], col = oi["blue"], lwd = 1.8)
  }
  if (z$main %in% c("stable node", "saddle")) {
    abline(h = 0, lwd = 2.5, col = oi["orange"]); abline(v = 0, lwd = 2.5, col = oi["orange"])
  }
  points(0, 0, pch = 21, bg = if (z$main == "saddle") "white" else "black", cex = 2)
}
dev.off()

## ---- fig04: trace-determinant plane -------------------------------------
Jsel <- function(x, a = 0.08) {
  u <- x[1]; v <- x[2]
  matrix(c(-1 + 2 * u * v, a + u^2, -2 * u * v, -(a + u^2)), 2, byrow = TRUE)
}
sel_fp <- c(0.6, 0.6 / (0.08 + 0.36))
svg("fig04-trace-det.svg", width = 8, height = 6)
par(mar = c(5, 5.2, 1, 1)); big()
plot(NA, xlim = c(-3, 3), ylim = c(-1, 2.5), xlab = expression("trace " * tau),
     ylab = expression("determinant " * Delta))
tt <- seq(-3, 3, length.out = 300)
polygon(c(-3, 3, 3, -3), c(-1, -1, 0, 0), col = "#F5E6D3", border = NA)
lines(tt, tt^2 / 4, lwd = 2.5, lty = 2)
abline(h = 0, lwd = 2); abline(v = 0, lwd = 2, col = oi["grey"])
text(0, -0.55, "saddle", cex = 1.4)
text(-2.5, 0.75, "stable\nnode", cex = 1.3)
text(2.5, 0.75, "unstable\nnode", cex = 1.3)
text(-1.1, 2.1, "stable\nspiral", cex = 1.3)
text(1.1, 2.1, "unstable\nspiral", cex = 1.3)
text(0.1, 1.4, "centre", cex = 1.2, srt = 90, adj = 0)
text(2.3, 1.6, expression(tau^2 == 4 * Delta), cex = 1.3)
pt <- function(A, col, lab, pos) {
  p <- c(sum(diag(A)), det(A)); points(p[1], p[2], pch = 21, bg = col, cex = 2.2, lwd = 1.5)
  text(p[1], p[2], lab, pos = pos, cex = 1.2)
}
pt(Jt(x_hi), "black", "toggle: U on / V on", 4)
pt(Jt(x_mid), "white", "toggle: middle", 4)
pt(Jsel(sel_fp), oi["orange"], "Sel'kov", 4)
dev.off()

## ---- fig05: the saddle up close ----------------------------------------
e <- eigen(Jt(x_mid))
svg("fig05-saddle.svg", width = 7, height = 6.5)
par(mar = c(5, 5.2, 1, 1)); big()
plot(NA, xlim = c(0.6, 1.8), ylim = c(0.6, 1.8), xlab = "u", ylab = "v", asp = 1)
field(ftog, c(0.65, 1.75), c(0.65, 1.75), nx = 12, len = 0.05)
abline(a = 0, b = 1, col = oi["orange"], lwd = 3)                         # stable (1, 1)
abline(a = 2 * x_mid[1], b = -1, col = oi["red"], lwd = 3, lty = 2)      # unstable (1, -1)
for (y0 in list(c(0.62, 0.70), c(0.70, 0.62), c(1.78, 1.70), c(1.70, 1.78),
                c(0.62, 0.64), c(0.64, 0.62), c(1.78, 1.76), c(1.76, 1.78)))
  traj(tog, y0, 12, oi["blue"], lwd = 2)
points(x_mid[1], x_mid[2], pch = 21, bg = "white", cex = 2.4, lwd = 2)
legend("bottomright", c(expression("stable: " * lambda == -2.19 * ", direction (1, 1)"),
                        expression("unstable: " * lambda == 0.19 * ", direction (1, -1)")),
       col = c(oi["orange"], oi["red"]), lwd = 3, lty = c(1, 2), cex = 1.05, bg = "white")
dev.off()

## ---- fig06: biased toggle: Newton basins vs ODE basins + separatrix -------
au <- 3; av <- 2.9
F2 <- function(x) tog(0, x, NULL, au, av)[[1]]
newt <- function(x) { for (k in 1:50) x <- x + solve(Jt(x, au, av), -F2(x)); x }
s_hi <- newt(c(3, 0.3)); s_mid <- newt(c(1.1, 1.3)); s_lo <- newt(c(0.4, 2.5))
g <- seq(0.01, 4, length.out = 300)
U <- matrix(rep(g, times = length(g)), length(g)); V <- matrix(rep(g, each = length(g)), length(g))
u <- U; v <- V
for (k in 1:60) {                         # vectorised Newton, Cramer's rule
  f1 <- au / (1 + v^2) - u; f2 <- av / (1 + u^2) - v
  b <- -2 * au * v / (1 + v^2)^2; c <- -2 * av * u / (1 + u^2)^2
  dt <- 1 - b * c
  du <- (f1 + b * f2) / dt; dv <- (c * f1 + f2) / dt
  u <- u + du; v <- v + dv
}
cls <- function(u, v) {
  out <- matrix(4, nrow(u), ncol(u)); ok <- is.finite(u) & is.finite(v)
  out[ok & abs(u - s_hi[1]) + abs(v - s_hi[2]) < 1e-6] <- 1
  out[ok & abs(u - s_mid[1]) + abs(v - s_mid[2]) < 1e-6] <- 2
  out[ok & abs(u - s_lo[1]) + abs(v - s_lo[2]) < 1e-6] <- 3
  out
}
Bn <- cls(u, v)
u <- U; v <- V
rk <- function(u, v) list(au / (1 + v^2) - u, av / (1 + u^2) - v)
h <- 0.05
for (k in 1:2000) {                        # vectorised RK4 to t = 100
  k1 <- rk(u, v); k2 <- rk(u + h / 2 * k1[[1]], v + h / 2 * k1[[2]])
  k3 <- rk(u + h / 2 * k2[[1]], v + h / 2 * k2[[2]]); k4 <- rk(u + h * k3[[1]], v + h * k3[[2]])
  u <- u + h / 6 * (k1[[1]] + 2 * k2[[1]] + 2 * k3[[1]] + k4[[1]])
  v <- v + h / 6 * (k1[[2]] + 2 * k2[[2]] + 2 * k3[[2]] + k4[[2]])
}
Bo <- ifelse(u > v, 1, 3)
es <- eigen(Jt(s_mid, au, av)); vs <- es$vectors[, es$values < 0]
sep <- rbind(ode(s_mid + 1e-6 * vs, seq(0, -7.5, by = -0.005), tog, NULL, au = au, av = av)[, 2:3],
             ode(s_mid - 1e-6 * vs, seq(0, -7.5, by = -0.005), tog, NULL, au = au, av = av)[, 2:3])
cols <- c(oi["green"], "grey85", oi["purple"], "black")
svg("fig06-separatrix.svg", width = 12, height = 5.8)
par(mfrow = c(1, 2), mar = c(5, 5.2, 3, 1)); big()
for (B in list(list(Bn, "Newton converges to ..."), list(Bo, "the cell (ODE) ends in ..."))) {
  image(g, g, B[[1]], col = cols, breaks = 0.5:4.5, xlab = "start u", ylab = "start v",
        main = B[[2]], asp = 1, useRaster = TRUE)
  points(rbind(s_hi, s_lo), pch = 21, bg = "black", cex = 2, lwd = 2)
  points(s_mid[1], s_mid[2], pch = 21, bg = "white", cex = 2, lwd = 2)
}
ord <- order(sep[, 1]); sep <- sep[ord, ]
sep <- sep[sep[, 1] >= 0 & sep[, 1] <= 4 & sep[, 2] >= 0 & sep[, 2] <= 4, ]
lines(sep, lwd = 3, col = oi["orange"])
segments(0, 0, 4, 4, lty = 3, lwd = 2)
legend(1.2, 4, c("U on", "middle", "V on", "separatrix", "u = v"), fill = c(cols[1:3], NA, NA),
       border = c("black", "black", "black", NA, NA), lty = c(NA, NA, NA, 1, 3),
       lwd = c(NA, NA, NA, 3, 2), col = c(NA, NA, NA, oi["orange"], "black"), cex = 1.1, bg = "white")
dev.off()
cat("biased toggle, Newton shares (U on, mid, V on, fail):", round(tabulate(Bn, 4) / length(Bn), 3), "\n")
cat("biased toggle, ODE shares (U on, V on):", round(c(mean(Bo == 1), mean(Bo == 3)), 3), "\n")

## ---- fig07: slow passage near the saddle ----------------------------------
svg("fig07-saddle-slow.svg", width = 9, height = 5.2)
par(mar = c(5, 5.2, 1, 1)); big()
tt <- seq(0, 100, by = 0.05)
plot(NA, xlim = c(0, 100), ylim = c(0, 3), xlab = "time (protein lifetimes)", ylab = "u (protein U)")
eps <- c(0, 1e-8, 1e-4, 1e-1)
lt <- c(1, 2, 4, 1)
cc <- c(oi["grey"], oi["blue"], oi["orange"], oi["red"])
for (i in seq_along(eps)) {
  o <- ode(c(2, 2 - eps[i]), tt, tog, NULL, rtol = 1e-10, atol = 1e-12)
  lines(o[, 1], o[, 2], lwd = 3, col = cc[i], lty = lt[i])
}
abline(h = x_mid[1], lty = 3)
text(98, x_mid[1] + 0.12, "middle state (saddle)", adj = 1, cex = 1.2)
legend("topleft", c("start (2, 2)", expression("start (2, 2 - " * 10^-8 * ")"),
                    expression("start (2, 2 - " * 10^-4 * ")"), "start (2, 1.9)"),
       col = cc, lty = lt, lwd = 3, cex = 1.15, bg = "white")
dev.off()

## ---- fig08: Sel'kov limit cycle -------------------------------------------
a <- 0.08; b <- 0.6
selkov <- function(t, x, p) list(c(-x[1] + a * x[2] + x[1]^2 * x[2], b - a * x[2] - x[1]^2 * x[2]))
fsel <- function(u, v) list(-u + a * v + u^2 * v, b - a * v - u^2 * v)
svg("fig08-selkov.svg", width = 13, height = 5.6)
par(mfrow = c(1, 2), mar = c(5, 5.2, 1, 1)); big()
plot(NA, xlim = c(0, 2.6), ylim = c(0, 3), xlab = "u (ADP)", ylab = "v (F6P)")
field(fsel, c(0.1, 2.5), c(0.1, 2.9), nx = 13, len = 0.035)
s <- seq(0, 2.6, length.out = 400)
lines(s, s / (a + s^2), col = oi["blue"], lwd = 3)
lines(s, b / (a + s^2), col = oi["red"], lwd = 3)
o1 <- traj(selkov, c(0.61, 1.36), 60, oi["orange"], lwd = 2, dt = 0.01)
o2 <- traj(selkov, c(2.5, 0.2), 60, oi["sky"], lwd = 2, dt = 0.01)
lc <- o1[o1[, 1] > 50, ]
lines(lc[, 2], lc[, 3], lwd = 4)
points(b, b / (a + b^2), pch = 21, bg = "white", cex = 2.2, lwd = 2)
legend("topright", c("du/dt = 0", "dv/dt = 0", "limit cycle"), col = c(oi["blue"], oi["red"], "black"),
       lwd = c(3, 3, 4), cex = 1.15, bg = "white")
plot(o1[, 1], o1[, 2], type = "l", lwd = 2.5, col = oi["orange"], xlab = "time",
     ylab = "u (ADP)", ylim = c(0, 2.6))
lines(o2[, 1], o2[, 2], lwd = 2.5, col = oi["sky"], lty = 2)
abline(h = b, lty = 3)
legend("topright", c("start (0.61, 1.36), near the fixed point", "start (2.5, 0.2)"),
       col = c(oi["orange"], oi["sky"]), lty = c(1, 2), lwd = 2.5, cex = 1.1, bg = "white")
dev.off()

## ---- fig09: bifurcation diagram in alpha ---------------------------------
al <- seq(1, 4, length.out = 400)
smid <- sapply(al, function(A) uniroot(function(s) s^3 + s - A, c(0, 3), tol = 1e-13)$root)
gmid <- 2 * smid^2 / (1 + smid^2)                    # off-diagonal size at u = v = s
ab <- al[al >= 2]
uhi <- (ab + sqrt(ab^2 - 4)) / 2; ulo <- (ab - sqrt(ab^2 - 4)) / 2
svg("fig09-bifurcation.svg", width = 13, height = 5.4)
par(mfrow = c(1, 2), mar = c(5, 5.2, 1, 1)); big()
plot(NA, xlim = c(1, 4), ylim = c(0, 3.8), xlab = expression("promoter strength " * alpha),
     ylab = "u at steady state")
lines(al[al <= 2], smid[al <= 2], lwd = 4)
lines(al[al >= 2], smid[al >= 2], lwd = 3, lty = 2)
lines(ab, uhi, lwd = 4, col = oi["green"]); lines(ab, ulo, lwd = 4, col = oi["purple"])
abline(v = 3, lty = 3, col = oi["grey"])
text(3.95, 3.55, "U on", adj = 1, cex = 1.3, col = oi["green"])
text(3.95, 0.05, "V on", adj = 1, cex = 1.3, col = oi["purple"])
text(3.95, 1.6, "middle (unstable)", adj = 1, cex = 1.3)
text(2, 0.5, expression(alpha == 2), cex = 1.3)
legend("topleft", c("stable", "unstable"), lty = c(1, 2), lwd = c(4, 3), cex = 1.2, bg = "white")
plot(NA, xlim = c(1, 4), ylim = c(-0.6, 0.4), xlab = expression("promoter strength " * alpha),
     ylab = expression("largest eigenvalue " * lambda[max]))
abline(h = 0, col = oi["grey"], lwd = 2)
lines(al[al <= 2], -1 + gmid[al <= 2], lwd = 4)
lines(al[al >= 2], -1 + gmid[al >= 2], lwd = 3, lty = 2)
lines(ab, -1 + 2 / ab, lwd = 4, col = oi["green"])
text(3.9, 0.3, "middle", adj = 1, cex = 1.3)
text(3.9, -0.25, "U on / V on: -1 + 2/α", adj = 1, cex = 1.3, col = oi["green"])
text(2.05, -0.55, "bifurcation at α = 2", adj = 0, cex = 1.3)
abline(v = 2, lty = 3)
dev.off()
