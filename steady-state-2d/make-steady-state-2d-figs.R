# Figures for the steady-state-2d/ unit ("Steady state in two dimensions").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd steady-state-2d && Rscript make-steady-state-2d-figs.R
#
# Base R graphics only (grDevices::svg) plus the Matrix package (ships with R).
# Colours: Okabe-Ito for lines, viridis (hcl.colors) for maps; both colour-blind safe.
# Running example: oxygen under a gas-impermeable lid. Chamber 4000 x 2000 um
# (top view), lid 100 um above the cells, open along all four edges to
# air-saturated medium (c0 = 200 uM). D = 2000 um^2/s (2e-9 m^2/s).
# Consumption per cell 50 amol/s (assumption); density rho in 1e5 cells/cm^2
# gives a volume rate R0 = 0.5 * rho uM/s. Michaelis-Menten Km = 1 uM.
# Grid h = 50 um: 79 x 39 interior points, x index i runs fastest (column-major).
# Makes:
#   fig01-chamber.svg      side and top view of the culture chamber under a lid
#   fig02-profiles-1d.svg  1D zero-order profiles across the chamber; c < 0 at high density
#   fig03-convergence.svg  1D linear uptake: error vs h, log-log, slope 2
#   fig04-stencil.svg      grid, numbering k = i + (j - 1) nx, five-point stencil
#   fig05-spy.svg          pattern of the 2D matrix (nx = 5, ny = 3) and its Kronecker blocks
#   fig06-ordering-bug.svg correct field vs the same solve with the Kronecker order swapped
#   fig07-oxygen-map.svg   zero-order oxygen maps for two densities
#   fig08-uptake.svg       zero-order, linear and Michaelis-Menten uptake vs c
#   fig09-newton.svg       Newton residuals and min(c): plain (wrong root) vs clipped
#   fig10-density.svg      Michaelis-Menten maps for three densities, hypoxic fraction vs density

library(Matrix)
oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999", black = "#000000")
um <- "µm"; uM <- "µM"
par_big <- function(...) par(cex.lab = 1.5, cex.axis = 1.35, cex.main = 1.5, ...)
pal <- hcl.colors(64, "viridis")

## ---- the running example (same code as on the slides) ----------------------

lap1 <- function(n) bandSparse(n, k = -1:1,
  diagonals = list(rep(-1, n - 1), rep(2, n), rep(-1, n - 1)))
D <- 2000; c0 <- 200; Km <- 1
setup <- function(h, Lx = 4000, Ly = 2000, swap = FALSE) {
  nx <- Lx / h - 1; ny <- Ly / h - 1
  A <- if (!swap) kronecker(Diagonal(ny), lap1(nx)) + kronecker(lap1(ny), Diagonal(nx))
       else       kronecker(lap1(nx), Diagonal(ny)) + kronecker(Diagonal(nx), lap1(ny))
  A <- forceSymmetric(D / h^2 * A)
  B <- matrix(0, nx, ny)
  B[1, ] <- B[1, ] + c0; B[nx, ] <- B[nx, ] + c0
  B[, 1] <- B[, 1] + c0; B[, ny] <- B[, ny] + c0
  list(A = A, b = D / h^2 * as.vector(B), nx = nx, ny = ny,
       x = h * (1:nx), y = h * (1:ny), h = h)
}
newton <- function(s, R0, clip = TRUE, start = c0, maxit = 50) {
  cv <- rep(start, s$nx * s$ny); res <- mn <- c()
  for (it in 1:maxit) {
    F <- as.vector(s$A %*% cv) + R0 * cv / (Km + cv) - s$b
    res <- c(res, sqrt(sum(F^2))); mn <- c(mn, min(cv))
    if (res[it] < 1e-8) break
    J <- s$A + Diagonal(x = R0 * Km / (Km + cv)^2)
    cv <- cv + as.vector(solve(J, -F))
    if (clip) cv <- pmax(cv, 0)
  }
  list(c = matrix(cv, s$nx, s$ny), res = res, min = mn)
}
s50 <- setup(50)
# map with the boundary values added, for plotting
full <- function(s, C) {
  list(x = c(0, s$x, s$x[s$nx] + s$h), y = c(0, s$y, s$y[s$ny] + s$h),
       z = rbind(c0, cbind(c0, C, c0), c0))
}
plotmap <- function(s, C, main, zlim = c(0, 200), neg = FALSE, hyp = TRUE) {
  f <- full(s, C)
  image(f$x, f$y, pmax(f$z, zlim[1]), zlim = zlim, col = pal, useRaster = TRUE,
        xlab = paste0("x (", um, ")"), ylab = paste0("y (", um, ")"), main = main)
  if (hyp) contour(f$x, f$y, f$z, levels = 10, add = TRUE, lwd = 3, col = oi["red"],
                   drawlabels = FALSE)
  if (neg && min(f$z) < 0) contour(f$x, f$y, f$z, levels = 0, add = TRUE, lwd = 3,
                                   lty = 2, col = "white", drawlabels = FALSE)
}
colbar <- function(zlim = c(0, 200), lab = paste0("O2 (", uM, ")")) {
  par(mar = c(4.8, 1, 3, 4.5))
  zz <- seq(zlim[1], zlim[2], length.out = 64)
  image(1, zz, matrix(zz, 1), col = pal, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1); mtext(lab, side = 4, line = 3.2, cex = 1.1)
}

## ---- fig01: the chamber ------------------------------------------------------

svg("fig01-chamber.svg", width = 12, height = 4.4)
par_big(mfrow = c(1, 2), mar = c(1, 1, 3, 1))
# side view
plot(NA, xlim = c(-0.6, 4.6), ylim = c(-0.6, 2.2), axes = FALSE, xlab = "", ylab = "",
     main = "side view (not to scale)")
rect(-0.6, -0.6, 4.6, 1.5, col = "#56B4E922", border = NA)
rect(-0.6, -0.6, 4.6, -0.3, col = "grey70", border = NA)
text(2, -0.45, "dish bottom", cex = 1.2)
rect(0, 0.75, 4, 0.95, col = "grey40", border = NA)
text(2, 1.12, "gas-impermeable lid", cex = 1.25)
set.seed(2)
xs <- seq(0.1, 3.9, by = 0.16)
points(xs, rep(-0.2, length(xs)), pch = 21, bg = oi["orange"], cex = 2.1)
text(2, 0.45, paste0("medium, h = 100 ", um), cex = 1.25)
arrows(-0.45, 1.3, -0.05, 0.2, length = 0.12, lwd = 3, col = oi["blue"])
arrows(4.45, 1.3, 4.05, 0.2, length = 0.12, lwd = 3, col = oi["blue"])
text(2, 1.85, paste0("air-saturated medium, c0 = 200 ", uM), col = oi["blue"], cex = 1.25)
text(2, 0.1, "cells consume O2", col = "#9A5B00", cex = 1.2)
# top view
s <- setup(50); res <- newton(s, 1)
f <- full(s, res$c)
par(mar = c(4.5, 4.8, 3, 1))
image(f$x, f$y, f$z, col = pal, useRaster = TRUE, zlim = c(0, 200),
      xlab = paste0("x (", um, ")"), ylab = paste0("y (", um, ")"),
      main = "top view: O2 under the lid")
contour(f$x, f$y, f$z, levels = 10, add = TRUE, lwd = 3, col = oi["red"], drawlabels = FALSE)
text(2000, 1000, "hypoxic\n(< 10 µM)", col = "white", cex = 1.3)
text(2000, 150, "edge: c = 200 µM", col = "black", cex = 1.2)
dev.off()

## ---- fig02: 1D profiles ------------------------------------------------------

svg("fig02-profiles-1d.svg", width = 9, height = 5.2)
par_big(mar = c(4.8, 5.2, 1, 1))
L <- 2000; yy <- seq(0, L, length.out = 401)
plot(NA, xlim = c(0, L), ylim = c(-60, 270), xlab = paste0("position across the chamber y (", um, ")"),
     ylab = paste0("O2 (", uM, ")"), las = 1)
rect(-100, -100, L + 100, 0, col = "grey90", border = NA)
abline(h = 0, col = "grey40")
text(L / 2, -35, "impossible: c < 0", cex = 1.3, col = "grey30")
R0s <- c(0.25, 0.5, 1); cols <- c(oi["sky"], oi["blue"], oi["red"])
for (j in 1:3) lines(yy, c0 - R0s[j] * yy * (L - yy) / (2 * D), lwd = 3, col = cols[j])
points(c(0, 500, 1000, 1500, 2000), c(200, 106.25, 75, 106.25, 200), pch = 21,
       bg = "white", col = oi["blue"], cex = 2, lwd = 2)
legend("top", ncol = 2, c(expression(0.5 %*% 10^5 ~ "cells/cm"^2), expression(1 %*% 10^5),
                     expression(2 %*% 10^5), "grid N = 4"),
       col = c(cols, oi["blue"]), lwd = c(3, 3, 3, NA), pch = c(NA, NA, NA, 21),
       pt.bg = "white", bty = "n", cex = 1.2)
dev.off()

## ---- fig03: 1D convergence ----------------------------------------------------

svg("fig03-convergence.svg", width = 8, height = 5.2)
par_big(mar = c(5.5, 7, 1, 1), mgp = c(4.5, 1, 0))
k <- 0.5 / 200; lam <- sqrt(D / k)
Ns <- c(10, 20, 40, 80, 160, 320); errs <- numeric(length(Ns))
for (m in seq_along(Ns)) {
  h <- L / Ns[m]; n <- Ns[m] - 1; y1 <- h * (1:n)
  A1 <- D / h^2 * lap1(n) + Diagonal(n, k)
  b1 <- numeric(n); b1[c(1, n)] <- D / h^2 * c0
  cn <- as.vector(solve(A1, b1))
  errs[m] <- max(abs(cn - c0 * cosh((y1 - L / 2) / lam) / cosh(L / (2 * lam))))
}
hs <- L / Ns
plot(hs, errs, log = "xy", pch = 19, cex = 1.6, col = oi["blue"], las = 1,
     xlab = paste0("grid spacing h (", um, ")"), ylab = paste0("max error (", uM, ")"))
lines(hs, errs, col = oi["blue"], lwd = 2)
lines(hs, 0.3 * errs[1] * (hs / hs[1])^2, lty = 2, lwd = 2, col = "grey40")
text(hs[3], errs[3] * 0.15, "slope 2", cex = 1.4, col = "grey30")
dev.off()

## ---- fig04: grid, numbering, stencil -----------------------------------------

svg("fig04-stencil.svg", width = 9.5, height = 5.4)
par_big(mar = c(3.5, 3.5, 0.5, 0.5))
nx <- 5; ny <- 3
plot(NA, xlim = c(-0.4, nx + 5.5), ylim = c(-0.4, ny + 1.4), axes = FALSE, asp = 1,
     xlab = "", ylab = "")
abline(h = 0:(ny + 1), col = "grey85"); abline(v = 0:(nx + 1), col = "grey85")
g <- expand.grid(i = 0:(nx + 1), j = 0:(ny + 1))
bd <- g$i == 0 | g$i == nx + 1 | g$j == 0 | g$j == ny + 1
points(g$i[bd], g$j[bd], pch = 22, bg = oi["blue"], cex = 2.4)
gi <- g[!bd, ]
kk <- gi$i + (gi$j - 1) * nx
st <- c(3 + nx * 1, 2 + nx * 1, 4 + nx * 1, 3, 3 + nx * 2)  # centre k = 8 and neighbours
bg <- ifelse(kk == 8, oi["red"], ifelse(kk %in% st, oi["orange"], "white"))
symbols(gi$i, gi$j, circles = rep(0.3, nrow(gi)), inches = FALSE, add = TRUE, bg = bg, lwd = 1.5)
text(gi$i, gi$j, kk, cex = 1.35)
mtext(paste0("x index i = 1 ... ", nx, "  (runs fastest)"), side = 1, line = 1.5, cex = 1.4)
mtext(paste0("y index j = 1 ... ", ny), side = 2, line = 1.5, cex = 1.4)
text(nx + 1.4, ny + 1.35, "boundary: known c", col = oi["blue"], cex = 1.2, pos = 2)
text(6.8, 2.0, "k = i + (j - 1) nx", cex = 1.4, pos = 4)
text(6.8, 1.5, "neighbours of k:", cex = 1.25, pos = 4)
text(6.8, 1.05, "k ± 1 (x),  k ± nx (y)", cex = 1.25, pos = 4)
dev.off()

## ---- fig05: the matrix pattern -------------------------------------------------

svg("fig05-spy.svg", width = 12, height = 4.6)
spy <- function(M, main) {
  M <- as.matrix(M); n <- nrow(M)
  plot(NA, xlim = c(0.5, n + 0.5), ylim = c(n + 0.5, 0.5), axes = FALSE, asp = 1,
       xlab = "", ylab = "", main = main)
  rect(0.5, n + 0.5, n + 0.5, 0.5, border = "grey60")
  w <- which(M != 0, arr.ind = TRUE)
  col <- ifelse(w[, 1] == w[, 2], oi["blue"], ifelse(abs(w[, 1] - w[, 2]) == 1, oi["orange"], oi["green"]))
  rect(w[, 2] - 0.45, w[, 1] - 0.45, w[, 2] + 0.45, w[, 1] + 0.45, col = col, border = NA)
}
par(mfrow = c(1, 3), mar = c(0.5, 0.5, 3, 0.5), cex.main = 1.6)
spy(kronecker(Diagonal(ny), lap1(nx)), "kronecker(I_ny, T_nx): x-neighbours")
abline(h = c(5.5, 10.5), v = c(5.5, 10.5), col = "grey60", lty = 3)
spy(kronecker(lap1(ny), Diagonal(nx)), "kronecker(T_ny, I_nx): y-neighbours")
abline(h = c(5.5, 10.5), v = c(5.5, 10.5), col = "grey60", lty = 3)
spy(kronecker(Diagonal(ny), lap1(nx)) + kronecker(lap1(ny), Diagonal(nx)), "sum: five-point Laplacian")
abline(h = c(5.5, 10.5), v = c(5.5, 10.5), col = "grey60", lty = 3)
dev.off()

## ---- fig06: the ordering bug ---------------------------------------------------

svg("fig06-ordering-bug.svg", width = 13, height = 4.3)
layout(matrix(1:3, 1), widths = c(1, 1, 0.22))
par_big(mar = c(4.8, 5, 3, 1))
sw <- setup(50, swap = TRUE)
Cok <- matrix(as.vector(solve(s50$A, s50$b - 0.5)), s50$nx, s50$ny)
Cbad <- matrix(as.vector(solve(sw$A, s50$b - 0.5)), s50$nx, s50$ny)
f1 <- full(s50, Cok); f2 <- full(s50, Cbad)
zl <- c(0, 300)
image(f1$x, f1$y, f1$z, zlim = zl, col = pal, useRaster = TRUE, main = "kron order matches as.vector",
      xlab = paste0("x (", um, ")"), ylab = paste0("y (", um, ")"))
image(f2$x, f2$y, pmin(f2$z, zl[2]), zlim = zl, col = pal, useRaster = TRUE, main = "kron order swapped (no error!)",
      xlab = paste0("x (", um, ")"), ylab = paste0("y (", um, ")"))
colbar(zl)
dev.off()

## ---- fig07: zero-order maps ------------------------------------------------------

svg("fig07-oxygen-map.svg", width = 13, height = 4.3)
layout(matrix(1:3, 1), widths = c(1, 1, 0.22))
par_big(mar = c(4.8, 5, 3, 1))
C1 <- matrix(as.vector(solve(s50$A, s50$b - 0.5)), s50$nx, s50$ny)
C2 <- matrix(as.vector(solve(s50$A, s50$b - 1)), s50$nx, s50$ny)
plotmap(s50, C1, sprintf("1e5 cells/cm2: min %.0f µM", min(C1)), hyp = FALSE)
plotmap(s50, C2, sprintf("2e5 cells/cm2: min %.1f µM", min(C2)), neg = TRUE, hyp = FALSE)
text(2000, 1000, "c < 0", col = "white", cex = 1.6)
colbar()
dev.off()

## ---- fig08: uptake models ----------------------------------------------------------

svg("fig08-uptake.svg", width = 12, height = 4.8)
par_big(mfrow = c(1, 2), mar = c(4.8, 5.5, 3, 1))
cc <- seq(0, 200, length.out = 801); R0 <- 1
plot(cc, rep(R0, length(cc)), type = "l", lwd = 3, col = oi["grey"], ylim = c(0, 1.3), las = 1,
     xlab = paste0("O2 (", uM, ")"), ylab = paste0("uptake r(c) (", uM, "/s)"),
     main = "full range")
lines(cc, R0 * cc / 200, lwd = 3, col = oi["green"], lty = 2)
lines(cc, R0 * cc / (Km + cc), lwd = 3, col = oi["blue"])
legend("bottomright", c("zero order: R0", "linear: (R0/c0) c", "Michaelis-Menten"),
       col = c(oi["grey"], oi["green"], oi["blue"]), lwd = 3, lty = c(1, 2, 1), bty = "n", cex = 1.2)
cs <- seq(0, 10, length.out = 401)
plot(cs, rep(R0, length(cs)), type = "l", lwd = 3, col = oi["grey"], ylim = c(0, 1.3), las = 1,
     xlab = paste0("O2 (", uM, ")"), ylab = paste0("uptake r(c) (", uM, "/s)"),
     main = "zoom: 0 to 10 µM")
lines(cs, R0 * cs / 200, lwd = 3, col = oi["green"], lty = 2)
lines(cs, R0 * cs / (Km + cs), lwd = 3, col = oi["blue"])
abline(v = Km, lty = 3); text(Km, 0.1, "Km", pos = 4, cex = 1.3)
dev.off()

## ---- fig09: Newton, plain vs clipped -------------------------------------------------

svg("fig09-newton.svg", width = 12, height = 4.8)
par_big(mfrow = c(1, 2), mar = c(5.5, 7, 3, 1), mgp = c(4.3, 1, 0))
np_ <- newton(s50, 1, clip = FALSE); nc_ <- newton(s50, 1, clip = TRUE)
plot(seq_along(np_$res) - 1, np_$res, log = "y", type = "b", pch = 19, lwd = 2, col = oi["red"],
     ylim = c(1e-12, 1e3), las = 1, xlab = "Newton iteration", ylab = "||F(c)|| (µM/s)",
     main = "residual")
lines(seq_along(nc_$res) - 1, nc_$res, type = "b", pch = 17, lwd = 2, col = oi["blue"])
legend("topright", c("plain Newton", "clipped at 0"), col = c(oi["red"], oi["blue"]),
       pch = c(19, 17), lwd = 2, bty = "n", cex = 1.2)
plot(seq_along(np_$min) - 1, np_$min, type = "b", pch = 19, lwd = 2, col = oi["red"],
     ylim = c(-40, 200), las = 1, xlab = "Newton iteration", ylab = paste0("min c (", uM, ")"),
     main = "smallest oxygen value")
lines(seq_along(nc_$min) - 1, nc_$min, type = "b", pch = 17, lwd = 2, col = oi["blue"])
abline(h = 0, lty = 3)
text(10, -13, sprintf("plain: %.1f µM", tail(np_$min, 1)), col = oi["red"], cex = 1.3)
text(7, 30, sprintf("clipped: %.2f µM", tail(nc_$min, 1)), col = oi["blue"], cex = 1.3)
dev.off()

## ---- fig10: density ---------------------------------------------------------------------

svg("fig10-density.svg", width = 13, height = 7.2)
layout(matrix(c(1, 2, 3, 4, 5, 5, 5, 6), 2, byrow = TRUE), widths = c(1, 1, 1, 0.25))
par_big(mar = c(4.5, 4.8, 3, 0.8))
for (rho in c(1, 2, 3)) {
  r <- newton(s50, 0.5 * rho)
  plotmap(s50, r$c, sprintf("%g e5 cells/cm2: %.0f %% hypoxic", rho, 100 * mean(r$c < 10)))
}
colbar()
rhos <- seq(0.25, 4, by = 0.25)
fr <- sapply(rhos, function(rho) mean(newton(s50, 0.5 * rho)$c < 10))
par(mar = c(4.8, 5.5, 1, 0.8))
plot(rhos, 100 * fr, type = "b", pch = 19, lwd = 2, col = oi["red"], las = 1,
     xlab = expression("cell density (" * 10^5 ~ "cells/cm"^2 * ")"),
     ylab = "hypoxic area (%)")
abline(h = 0, lty = 3)
plot.new()
dev.off()

cat("done\n")
