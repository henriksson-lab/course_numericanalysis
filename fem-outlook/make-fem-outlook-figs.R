# Figures for the fem-outlook/ unit ("Beyond rectangles: finite elements and flow").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd fem-outlook && Rscript make-fem-outlook-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito (colour-blind
# safe). Runs in a few seconds. Units: um and s throughout.
#
# Running examples (our own round numbers; all "data" are synthetic):
#   * culture chamber 600 x 400 um with a 3 x 2 array of round pillars,
#     radius 50 um, centres 200 um apart.
#   * drug diffusing into a hydrogel slab of thickness L = 1000 um with uptake
#     by embedded cells: -(D u')' + k u = 0, D = 100 um^2/s, k = 0.04 /s,
#     u(0) = 1 (bath), u'(L) = 0 (dish bottom). Length scale lambda = 50 um.
#   * Poiseuille flow between plates 100 um apart, mean speed 1000 um/s
#     (the T-sensor numbers of "Flow plus diffusion").
# Makes:
#   fig01-staircase.svg   chamber with pillars on a square grid: staircase walls
#   fig02-finite-volume.svg  conservation per cell: 1D row and a 2D triangle
#   fig03-hats.svg        hat functions on a non-uniform mesh, and their sum
#   fig04-assembly.svg    2 x 2 element matrices adding up to a tridiagonal one
#   fig05-slab.svg        hydrogel slab: exact profile, uniform vs graded mesh
#   fig06-error.svg       max nodal error vs number of elements, log-log
#   fig07-negative.svg    too-coarse mesh: negative concentration
#   fig08-mesh2d.svg      triangle mesh of the pillar chamber, graded at pillars
#   fig09-poiseuille.svg  Poiseuille profile across a channel, FEM nodes on it

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999", black = "#000000")
um <- "µm"
par_big <- function(...) par(cex.lab = 1.5, cex.axis = 1.35, cex.main = 1.5, ...)

## ---- solver (same code as on the slides) -------------------------------------

thomas <- function(a, b, c, d) {                 # from "Implicit schemes ..."
  n <- length(d)
  for (i in 2:n) { w <- a[i] / b[i - 1]; b[i] <- b[i] - w * c[i - 1]; d[i] <- d[i] - w * d[i - 1] }
  x <- numeric(n); x[n] <- d[n] / b[n]
  for (i in (n - 1):1) x[i] <- (d[i] - c[i] * x[i + 1]) / b[i]
  x
}
assemble <- function(x, D, k, f = 0, lumped = FALSE) {
  n <- length(x); h <- diff(x)
  D <- rep_len(D, n - 1); k <- rep_len(k, n - 1); f <- rep_len(f, n - 1)
  di <- numeric(n); off <- numeric(n - 1); b <- numeric(n)
  for (e in 1:(n - 1)) {
    kd <- D[e] / h[e] + k[e] * h[e] * (if (lumped) 1 / 2 else 1 / 3)
    ko <- -D[e] / h[e] + k[e] * h[e] * (if (lumped) 0 else 1 / 6)
    di[e:(e + 1)] <- di[e:(e + 1)] + kd
    off[e] <- off[e] + ko
    b[e:(e + 1)] <- b[e:(e + 1)] + f[e] * h[e] / 2
  }
  list(di = di, off = off, b = b)
}
fem1d <- function(x, D, k, f = 0, u0 = 1, lumped = FALSE) {
  S <- assemble(x, D, k, f, lumped); off <- S$off
  rhs <- S$b[-1]; rhs[1] <- rhs[1] - off[1] * u0
  c(u0, thomas(off, S$di[-1], c(off[-1], 0), rhs))
}

L <- 1000; D <- 100; k <- 0.04; lam <- sqrt(D / k)
exact <- function(x) cosh((L - x) / lam) / cosh(L / lam)

## ---- chamber geometry ----------------------------------------------------------

W <- 600; H <- 400; R <- 50
cx <- c(100, 300, 500); cy <- c(100, 300)
centres <- expand.grid(x = cx, y = cy)
in_pillar <- function(x, y) {
  d2 <- outer(x, centres$x, "-")^2 + outer(y, centres$y, "-")^2
  apply(d2 <= R^2, 1, any)
}
draw_pillars <- function(col = oi["grey"], border = oi["black"]) {
  th <- seq(0, 2 * pi, length.out = 200)
  for (i in seq_len(nrow(centres)))
    polygon(centres$x[i] + R * cos(th), centres$y[i] + R * sin(th),
            col = col, border = border, lwd = 2)
}

## ---- fig01: staircase ------------------------------------------------------------

svg("fig01-staircase.svg", width = 12, height = 4.6)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
plot(NA, xlim = c(0, W), ylim = c(0, H), asp = 1, axes = FALSE, xlab = "", ylab = "")
title("The chamber: round pillars", cex.main = 1.6)
rect(0, 0, W, H, col = "#EAF4FB", border = oi["black"], lwd = 2)
draw_pillars()
arrows(-5, 200, 60, 200, length = 0.12, lwd = 3, col = oi["blue"])
text(300, -12, paste("600 ×", 400, um, ", pillars r = 50", um), cex = 1.3, xpd = NA)

hgrid <- 25
plot(NA, xlim = c(0, W), ylim = c(0, H), asp = 1, axes = FALSE, xlab = "", ylab = "")
title(paste0("Square grid, Δx = ", hgrid, " ", um, ": staircase walls"), cex.main = 1.6)
xs <- seq(hgrid / 2, W - hgrid / 2, by = hgrid); ys <- seq(hgrid / 2, H - hgrid / 2, by = hgrid)
G <- expand.grid(x = xs, y = ys)
solid <- in_pillar(G$x, G$y)
rect(G$x - hgrid / 2, G$y - hgrid / 2, G$x + hgrid / 2, G$y + hgrid / 2,
     col = ifelse(solid, oi["grey"], "#EAF4FB"), border = "white", lwd = 0.8)
th <- seq(0, 2 * pi, length.out = 200)
for (i in seq_len(nrow(centres)))
  lines(centres$x[i] + R * cos(th), centres$y[i] + R * sin(th), col = oi["red"], lwd = 2.5)
rect(0, 0, W, H, border = oi["black"], lwd = 2)
dev.off()

## ---- fig02: finite volumes -----------------------------------------------------

svg("fig02-finite-volume.svg", width = 12, height = 4.2)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
plot(NA, xlim = c(0, 5), ylim = c(0, 2.2), axes = FALSE, xlab = "", ylab = "")
title("1D: cells of any width", cex.main = 1.6)
edges <- c(0.2, 1.0, 2.6, 3.4, 4.8)
for (i in 1:4) rect(edges[i], 0.6, edges[i + 1], 1.6,
                    col = if (i == 2) "#FCE9C8" else "#EAF4FB", border = oi["black"], lwd = 2)
arrows(0.75, 1.1, 1.3, 1.1, lwd = 3, length = 0.12, col = oi["blue"])
arrows(2.3, 1.1, 2.9, 1.1, lwd = 3, length = 0.12, col = oi["blue"])
arrows(1.8, 1.4, 1.8, 2.0, lwd = 3, length = 0.12, col = oi["red"])
text(1.0, 0.4, expression(J["in"]), cex = 1.5, col = oi["blue"])
text(2.6, 0.4, expression(J["out"]), cex = 1.5, col = oi["blue"])
text(2.0, 2.1, "uptake", cex = 1.4, col = oi["red"], adj = 0)
text(1.8, 0.8, "cell i", cex = 1.3)

plot(NA, xlim = c(0, 4), ylim = c(0, 3), asp = 1, axes = FALSE, xlab = "", ylab = "")
title("2D: the same balance on a triangle", cex.main = 1.6)
P <- rbind(c(0.8, 0.4), c(3.4, 0.7), c(1.8, 2.6))
polygon(P, col = "#FCE9C8", border = oi["black"], lwd = 2)
for (s in 1:3) {
  a <- P[s, ]; b <- P[s %% 3 + 1, ]; m <- (a + b) / 2
  t <- b - a; nrm <- c(t[2], -t[1]) / sqrt(sum(t^2))
  arrows(m[1] - 0.35 * nrm[1], m[2] - 0.35 * nrm[2], m[1] + 0.35 * nrm[1], m[2] + 0.35 * nrm[2],
         lwd = 3, length = 0.12, col = oi["blue"])
}
text(2.0, 1.25, "sum of 3 face fluxes\n= uptake inside", cex = 1.25)
dev.off()

## ---- fig03: hat functions --------------------------------------------------------

svg("fig03-hats.svg", width = 12, height = 4.6)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1)); par_big()
xn <- c(0, 1, 1.5, 2, 3, 5)
cols <- c(oi["blue"], oi["orange"], oi["green"], oi["red"], oi["purple"], oi["sky"])
plot(NA, xlim = range(xn), ylim = c(0, 1.15), xlab = "x", ylab = expression(phi[j](x)))
title("Hat functions on a non-uniform mesh", cex.main = 1.5)
for (j in seq_along(xn)) {
  y <- as.numeric(seq_along(xn) == j)
  lines(xn, y, col = cols[j], lwd = 3)
}
points(xn, rep(0, length(xn)), pch = 19, cex = 1.4)
text(1.5, 1.08, expression(phi[2]), cex = 1.4, col = oi["green"])
uj <- c(1, 0.45, 0.3, 0.2, 0.08, 0.02)
xf <- seq(0, 5, length.out = 400)
plot(xf, exp(-0.8 * xf), type = "l", lwd = 2, lty = 2, col = oi["grey"],
     xlab = "x", ylab = "u", ylim = c(0, 1.15))
title(expression(u[h](x) == sum(u[j] * phi[j](x), j)), cex.main = 1.5)
for (j in seq_along(xn)) {
  y <- uj[j] * as.numeric(seq_along(xn) == j)
  lines(xn, y, col = cols[j], lwd = 1.5)
}
lines(xn, uj, lwd = 3.5, col = oi["black"])
points(xn, uj, pch = 19, cex = 1.4)
legend("topright", c("exact u", expression(u[h]~"(piecewise linear)")), lty = c(2, 1),
       lwd = c(2, 3.5), col = c(oi["grey"], oi["black"]), bty = "n", cex = 1.3)
dev.off()

## ---- fig04: assembly -----------------------------------------------------------

svg("fig04-assembly.svg", width = 11, height = 4.8)
par(mar = c(1, 1, 3, 1))
n <- 6
plot(NA, xlim = c(0, 15.5), ylim = c(-0.5, n + 0.5), axes = FALSE, xlab = "", ylab = "", asp = 1)
title("Each element adds a 2 × 2 block; neighbours overlap on the diagonal", cex.main = 1.5)
ecol <- c(oi["blue"], oi["orange"], oi["green"], oi["red"], oi["purple"])
## element blocks, offset a little so overlap is visible
for (e in 1:(n - 1)) {
  x0 <- e - 1; y0 <- n - e - 1
  rect(x0 + 0.06, y0 + 0.06, x0 + 1.94, y0 + 1.94, border = ecol[e], lwd = 3.5)
}
for (i in 1:n) for (j in 1:n) if (abs(i - j) <= 1)
  text(j - 0.5, n - i + 0.5, if (i == j) "•" else "•", cex = 2)
rect(0, 0, n, n, border = oi["black"], lwd = 1.5)
text(n / 2, -0.4, "assembled matrix: tridiagonal", cex = 1.3)
## element matrix
text(11, 5.2, "element e, length h_e:", cex = 1.35)
text(11, 3.9, expression(frac(D, h[e]) * bgroup("(", atop(" 1  -1", "-1   1"), ")")), cex = 1.6)
text(11, 2.3, expression(+ frac(k * h[e], 6) * bgroup("(", atop("2  1", "1  2"), ")")), cex = 1.6)
text(11, 0.8, "added into rows/columns e, e+1", cex = 1.3)
dev.off()

## ---- fig05: slab profiles --------------------------------------------------------

xu <- seq(0, L, length.out = 11); xg <- L * (0:10 / 10)^2
uu <- fem1d(xu, D, k); ug <- fem1d(xg, D, k)
svg("fig05-slab.svg", width = 12, height = 4.8)
par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3, 1)); par_big()
xf <- seq(0, L, length.out = 1000)
for (zoom in c(FALSE, TRUE)) {
  xl <- if (zoom) c(0, 300) else c(0, L)
  plot(xf, exact(xf), type = "l", lwd = 2.5, col = oi["black"], xlim = xl, ylim = c(-0.02, 1),
       xlab = paste0("depth x (", um, ")"), ylab = "drug u / bath")
  title(if (zoom) "Top 300 µm" else "Whole slab, 10 elements each", cex.main = 1.5)
  lines(xu, uu, col = oi["orange"], lwd = 2.5, type = "o", pch = 19, cex = 1.2)
  lines(xg, ug, col = oi["blue"], lwd = 2.5, type = "o", pch = 17, cex = 1.2)
  if (!zoom) legend("topright", c("exact", "uniform mesh", "graded mesh"),
                    col = c(oi["black"], oi["orange"], oi["blue"]), lwd = 2.5,
                    pch = c(NA, 19, 17), bty = "n", cex = 1.3)
}
dev.off()

## ---- fig06: error vs N -------------------------------------------------------------

Ns <- c(5, 10, 20, 40, 80, 160, 320)
err <- t(sapply(Ns, function(N) {
  xu <- seq(0, L, length.out = N + 1); xg <- L * (0:N / N)^2
  c(max(abs(fem1d(xu, D, k) - exact(xu))), max(abs(fem1d(xg, D, k) - exact(xg))))
}))
svg("fig06-error.svg", width = 8, height = 5.2)
par(mar = c(4.8, 5, 2, 1)); par_big()
plot(Ns, err[, 1], log = "xy", type = "o", pch = 19, lwd = 2.5, col = oi["orange"], cex = 1.3,
     xlab = "number of elements N", ylab = "max nodal error", ylim = range(err))
lines(Ns, err[, 2], type = "o", pch = 17, lwd = 2.5, col = oi["blue"], cex = 1.3)
lines(Ns, 0.5 * err[3, 1] * (Ns / Ns[3])^-2, lty = 2, col = oi["grey"], lwd = 2)
text(160, 0.5 * err[3, 1] * (160 / Ns[3])^-2, "slope -2", pos = 1, cex = 1.3, col = oi["grey"])
legend("topright", c("uniform", "graded  x = L (i/N)²"), col = c(oi["orange"], oi["blue"]),
       pch = c(19, 17), lwd = 2.5, bty = "n", cex = 1.3)
dev.off()

## ---- fig07: negative concentration --------------------------------------------------

svg("fig07-negative.svg", width = 12, height = 4.6)
par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3, 1)); par_big()
xf <- seq(0, 600, length.out = 600)
x4 <- seq(0, L, length.out = 5); x8 <- seq(0, L, length.out = 9)
plot(xf, exact(xf), type = "l", lwd = 2.5, xlim = c(0, 600), ylim = c(-0.12, 1),
     xlab = paste0("depth x (", um, ")"), ylab = "drug u / bath")
abline(h = 0, col = oi["grey"])
title("Elements too long: h = 250 µm > 2.45 λ", cex.main = 1.5)
lines(x4, fem1d(x4, D, k), type = "o", pch = 19, lwd = 2.5, col = oi["red"], cex = 1.3)
lines(x4, fem1d(x4, D, k, lumped = TRUE), type = "o", pch = 15, lwd = 2.5, col = oi["green"], cex = 1.3)
legend("topright", c("exact", "FEM, 4 elements", "lumped (= FD)"),
       col = c(oi["black"], oi["red"], oi["green"]), lwd = 2.5, pch = c(NA, 19, 15), bty = "n", cex = 1.3)
hs <- seq(10, 400, by = 2)
plot(hs / lam, -D / hs + k * hs / 6, type = "l", lwd = 3, col = oi["blue"],
     xlab = "h / λ", ylab = "off-diagonal entry")
abline(h = 0, col = oi["grey"]); abline(v = sqrt(6), lty = 2, lwd = 2, col = oi["red"])
text(sqrt(6), 1.2, expression(sqrt(6)), pos = 4, cex = 1.4, col = oi["red"])
title("Sign of  -D/h + kh/6", cex.main = 1.5)
dev.off()

## ---- fig08: 2D triangle mesh -----------------------------------------------------

## O-grid around each pillar: blend circle and square cell boundary, graded
## towards the pillar; split each quadrilateral into two triangles.
svg("fig08-mesh2d.svg", width = 12, height = 4.6)
par(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
mesh_cell <- function(x0, y0, a = 100, r = R, M = 32, Rn = 6, grade = 1.6) {
  th <- 2 * pi * (0:(M - 1)) / M + pi / M
  s <- ((0:Rn) / Rn)^grade
  sq <- cbind(cos(th), sin(th)) / pmax(abs(cos(th)), abs(sin(th))) * a
  ci <- cbind(cos(th), sin(th)) * r
  X <- outer(1 - s, ci[, 1]) + outer(s, sq[, 1]) + x0
  Y <- outer(1 - s, ci[, 2]) + outer(s, sq[, 2]) + y0
  tris <- list()
  for (i in 1:Rn) for (j in 1:M) {
    j2 <- j %% M + 1
    tris[[length(tris) + 1]] <- cbind(c(X[i, j], X[i + 1, j], X[i + 1, j2]),
                                      c(Y[i, j], Y[i + 1, j], Y[i + 1, j2]))
    tris[[length(tris) + 1]] <- cbind(c(X[i, j], X[i + 1, j2], X[i, j2]),
                                      c(Y[i, j], Y[i + 1, j2], Y[i, j2]))
  }
  tris
}
alltris <- list()
for (i in seq_len(nrow(centres))) alltris <- c(alltris, mesh_cell(centres$x[i], centres$y[i]))
plot(NA, xlim = c(0, W), ylim = c(0, H), asp = 1, axes = FALSE, xlab = "", ylab = "")
title(paste(length(alltris), "triangles, fitted to the pillars"), cex.main = 1.6)
for (tr in alltris) polygon(tr, col = "#EAF4FB", border = oi["blue"], lwd = 0.6)
draw_pillars()
rect(0, 0, W, H, border = oi["black"], lwd = 2)
## zoom on one pillar
plot(NA, xlim = c(200, 400), ylim = c(200, 400), asp = 1, axes = FALSE, xlab = "", ylab = "")
title("Zoom: small triangles at the wall", cex.main = 1.6)
for (tr in mesh_cell(300, 300)) polygon(tr, col = "#EAF4FB", border = oi["blue"], lwd = 1.2)
th <- seq(0, 2 * pi, length.out = 200)
polygon(300 + R * cos(th), 300 + R * sin(th), col = oi["grey"], border = oi["black"], lwd = 2)
dev.off()
cat("triangles:", length(alltris), " grid cells in pillars (fig01):", sum(solid), "\n")

## ---- fig09: Poiseuille ---------------------------------------------------------------

w <- 100; ubar <- 1000; G <- 12 * ubar / w^2           # G = pressure gradient / viscosity
yh <- seq(0, w / 2, length.out = 6)
uh <- fem1d(yh, D = 1, k = 0, f = G, u0 = 0)            # half channel, symmetry at centre
svg("fig09-poiseuille.svg", width = 11, height = 4.6)
par(mfrow = c(1, 2), mar = c(4.5, 4.8, 3, 1)); par_big()
yf <- seq(0, w, length.out = 300)
plot(G / 2 * yf * (w - yf), yf, type = "l", lwd = 3, col = oi["blue"],
     xlab = paste0("speed (", um, "/s)"), ylab = paste0("y across channel (", um, ")"),
     xlim = c(0, 1600))
points(uh, yh, pch = 19, cex = 1.5, col = oi["red"])
abline(v = ubar, lty = 2, col = oi["grey"], lwd = 2)
text(ubar, 95, "mean", pos = 2, cex = 1.3, col = oi["grey"])
title("Straight channel: a parabola", cex.main = 1.5)
legend("bottomright", c("exact", "FEM nodes (5 elements)"), col = c(oi["blue"], oi["red"]),
       lwd = c(3, NA), pch = c(NA, 19), bty = "n", cex = 1.2)
## sketch: pillars, no formula
par(mar = c(1, 1, 3, 1))
plot(NA, xlim = c(0, W), ylim = c(0, H), asp = 1, axes = FALSE, xlab = "", ylab = "")
title("Around pillars: no formula", cex.main = 1.5)
rect(0, 0, W, H, col = "#EAF4FB", border = oi["black"], lwd = 2)
draw_pillars()
for (yy in c(20, 200, 380)) arrows(10, yy, 60, yy, length = 0.1, lwd = 2.5, col = oi["blue"])
text(300, 200, "?", cex = 4, col = oi["red"])
dev.off()

cat("uniform/graded max err N=10:", err[2, ], "\n")
print(round(uh, 6))
