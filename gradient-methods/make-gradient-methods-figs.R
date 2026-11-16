# Figures for the gradient-methods/ unit ("Gradient descent, Newton and
# quasi-Newton").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd gradient-methods && Rscript make-gradient-methods-figs.R
#
# Base R graphics only (grDevices::svg) plus numDeriv. Colours are Okabe-Ito
# (colour-blind safe). Running example: a 4-parameter logistic (4PL)
# dose-response curve fitted to 18 simulated wells (9 doses x 2), parameters
# p = (bottom, top, EC50 [uM], hill). Makes:
#   fig01-dose-response.svg   data and fitted 4PL curve on a log dose axis
#   fig02-gradient.svg        SSE contours in (EC50, hill) with -gradient arrows
#   fig03-zigzag.svg          steepest descent on an elongated valley; error vs kappa
#   fig04-convergence.svg     4PL: steepest descent vs Newton vs BFGS
#   fig05-nelder-mead.svg     Nelder-Mead moves and the points optim evaluates
#   fig06-units-plateau.svg   SSE vs EC50 with dose in nM: the start on a plateau

library(numDeriv)

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")
open_svg <- function(file, w, h) {
  svg(file, width = w, height = h)
  par(mar = c(4.6, 4.8, 1.5, 1), cex.lab = 1.5, cex.axis = 1.3, las = 1)
}

# ---- data and model (same numbers as on the slides) ------------------------
# simulated with set.seed(27), p = (5, 100, 1.5, 1.2), noise sd 4, rounded
dose <- rep(10^seq(-2, 2, by = 0.5), each = 2)
y <- c(12.9, 9.8, 2.9, 0.1, 4.2, 9.7, 17.7, 22.3, 49.7, 42.1,
       67.3, 72.6, 97.4, 91.8, 94.6, 93.3, 92.9, 95.1)
f4 <- function(p, c) p[1] + (p[2] - p[1]) / (1 + (p[3] / c)^p[4])
sse <- function(p) sum((y - f4(p, dose))^2)
grad_sse <- function(p) {
  u <- (p[3] / dose)^p[4]; D <- 1 + u
  r <- y - f4(p, dose)
  J <- cbind(u / D, 1 / D,
             -(p[2] - p[1]) * p[4] * u / (p[3] * D^2),
             -(p[2] - p[1]) * u * log(p[3] / dose) / D^2)
  -2 * drop(t(J) %*% r)
}
hess <- function(p) jacobian(grad_sse, p)
p0 <- c(0, 100, 1, 1)
pstar <- c(5.393180177899748, 96.173705843036529, 1.270980009586612,
           1.239553930417529)
nrm <- function(v) sqrt(sum(v^2))

# ---- fig01: data and fit ----------------------------------------------------
open_svg("fig01-dose-response.svg", 8, 5)
cc <- 10^seq(-2.3, 2.3, length.out = 300)
plot(dose, y, log = "x", pch = 19, cex = 1.5, col = oi["blue"],
     xlab = expression("dose (" * mu * "M)"), ylab = "response (% of max)",
     ylim = c(-5, 110), xaxt = "n")
axis(1, at = 10^(-2:2), labels = c("0.01", "0.1", "1", "10", "100"))
lines(cc, f4(pstar, cc), lwd = 3, col = oi["orange"])
abline(h = pstar[1:2], lty = 2, col = oi["grey"], lwd = 1.5)
mid <- (pstar[1] + pstar[2]) / 2
segments(pstar[3], -10, pstar[3], mid, lty = 3, lwd = 2)
segments(1e-3, mid, pstar[3], mid, lty = 3, lwd = 2)
text(0.006, pstar[2] + 6, "top", cex = 1.3, adj = 0)
text(0.006, pstar[1] - 6, "bottom", cex = 1.3, adj = 0)
text(pstar[3] * 1.2, 2, "EC50", cex = 1.3, adj = 0)
text(5, 40, "slope set\nby hill", cex = 1.3, adj = 0)
dev.off()

# ---- fig02: gradient arrows on a 2-D slice ---------------------------------
open_svg("fig02-gradient.svg", 7.5, 5.6)
ec <- seq(0.6, 2.0, length.out = 120); hl <- seq(0.7, 1.9, length.out = 120)
Z <- outer(ec, hl, Vectorize(function(a, b) sse(c(pstar[1:2], a, b))))
contour(ec, hl, Z, levels = c(320, 350, 400, 500, 700, 1000, 1500, 2000, 3000),
        xlab = expression("EC50 (" * mu * "M)"), ylab = "hill", asp = 1,
        labcex = 1, col = oi["grey"], lwd = 1.5, xlim = c(0.6, 2.0), ylim = c(0.7, 1.9))
pts <- as.matrix(expand.grid(seq(0.75, 1.85, by = 0.275), seq(0.8, 1.8, by = 0.25)))
for (i in seq_len(nrow(pts))) {
  g <- grad_sse(c(pstar[1:2], pts[i, ]))[3:4]
  d <- -g / nrm(g) * 0.12
  arrows(pts[i, 1], pts[i, 2], pts[i, 1] + d[1], pts[i, 2] + d[2],
         length = 0.08, lwd = 2.2, col = oi["blue"])
}
points(pstar[3], pstar[4], pch = 4, cex = 2, lwd = 3, col = oi["red"])
dev.off()

# ---- fig03: zig-zag on an elongated valley ---------------------------------
sd_exact <- function(kap, x, n) {      # f = (x1^2 + kap x2^2)/2, exact line search
  P <- matrix(x, nrow = 1)
  for (k in 1:n) {
    g <- c(x[1], kap * x[2])
    a <- sum(g^2) / (g[1]^2 + kap * g[2]^2)
    x <- x - a * g; P <- rbind(P, x)
  }
  P
}
svg("fig03-zigzag.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 2.2, 1), cex.lab = 1.5, cex.axis = 1.3,
    las = 1, cex.main = 1.4)
kap <- 10
xs <- seq(-11, 11, length.out = 200); ys <- seq(-3, 3, length.out = 200)
contour(xs, ys, outer(xs, ys, function(a, b) (a^2 + kap * b^2) / 2),
        levels = c(1, 5, 15, 30, 50), drawlabels = FALSE, asp = 1,
        col = oi["grey"], lwd = 1.5, xlab = expression(x[1]), ylab = expression(x[2]),
        main = expression("steepest descent, " * kappa * " = 10"))
P <- sd_exact(kap, c(10, 1), 12)
lines(P, type = "o", pch = 19, col = oi["red"], lwd = 2.5, cex = 0.9)
points(0, 0, pch = 4, cex = 2, lwd = 3)
plot(NA, xlim = c(0, 100), ylim = c(1e-12, 10), log = "y", yaxt = "n",
     xlab = "iteration", ylab = expression("error " * group("||", x[k] - x^"*", "||")),
     main = "exact line search, start (1, 1)")
axis(2, at = 10^c(-12, -8, -4, 0), labels = expression(10^-12, 10^-8, 10^-4, 10^0))
cols <- c(oi["blue"], oi["orange"], oi["red"])
for (i in 1:3) {
  kk <- c(10, 100, 1000)[i]
  P <- sd_exact(kk, c(1, 1), 100)
  lines(0:100, apply(P, 1, nrm), lwd = 3, col = cols[i])
}
legend("bottomleft", legend = expression(kappa == 10, kappa == 100, kappa == 1000),
       col = cols, lwd = 3, cex = 1.3, bty = "n")
dev.off()

# ---- fig04: convergence on the 4PL -----------------------------------------
sd_path <- function(x, tol = 1e-4, maxit = 5000) {
  P <- matrix(x, nrow = 1)
  for (k in 1:maxit) {
    g <- grad_sse(x); if (nrm(g) < tol) break
    a <- 1
    while (!isTRUE(sse(x - a * g) <= sse(x) - 1e-4 * a * sum(g^2))) a <- a / 2
    x <- x - a * g; P <- rbind(P, x)
  }
  P
}
newton_path <- function(x, tol = 1e-4, maxit = 50) {
  P <- matrix(x, nrow = 1)
  for (k in 1:maxit) {
    g <- grad_sse(x); if (nrm(g) < tol) break
    p <- solve(hess(x), -g)
    a <- 1
    while (!isTRUE(sse(x + a * p) <= sse(x) + 1e-4 * a * sum(g * p))) a <- a / 2
    x <- x + a * p; P <- rbind(P, x)
  }
  P
}
bfgs_path <- function(x, tol = 1e-4, maxit = 200) {
  n <- length(x); B <- diag(n); g <- grad_sse(x); P <- matrix(x, nrow = 1)
  for (k in 1:maxit) {
    p <- -drop(B %*% g)
    a <- 1
    while (!isTRUE(sse(x + a * p) <= sse(x) + 1e-4 * a * sum(g * p))) a <- a / 2
    s <- a * p; x <- x + s; gn <- grad_sse(x); yv <- gn - g; g <- gn
    P <- rbind(P, x)
    if (nrm(g) < tol) break
    if (sum(yv * s) > 0) {
      rho <- 1 / sum(yv * s); I <- diag(n)
      B <- (I - rho * s %*% t(yv)) %*% B %*% (I - rho * yv %*% t(s)) + rho * s %*% t(s)
    }
  }
  P
}
Psd <- sd_path(p0); Pn <- newton_path(p0); Pb <- bfgs_path(p0)
err <- function(P) pmax(apply(sweep(P, 2, pstar), 1, nrm), 1e-14)
cat("iterations: SD", nrow(Psd) - 1, " Newton", nrow(Pn) - 1, " BFGS", nrow(Pb) - 1, "\n")
svg("fig04-convergence.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 2.2, 1), cex.lab = 1.5, cex.axis = 1.3,
    las = 1, cex.main = 1.4)
yl <- c(1e-10, 100)
plot(NA, xlim = c(0, 30), ylim = yl, log = "y", yaxt = "n", xlab = "iteration",
     ylab = expression("error " * group("||", p[k] - p^"*", "||")), main = "first 30 iterations")
axis(2, at = 10^c(-10, -6, -2, 2), labels = expression(10^-10, 10^-6, 10^-2, 10^2))
lines(0:(nrow(Psd) - 1), err(Psd), lwd = 3, col = oi["red"])
lines(0:(nrow(Pn) - 1), err(Pn), type = "o", pch = 19, lwd = 3, col = oi["blue"])
lines(0:(nrow(Pb) - 1), err(Pb), type = "o", pch = 17, lwd = 3, col = oi["orange"])
legend("bottomleft", c("steepest descent", "Newton", "BFGS"), lwd = 3, pch = c(NA, 19, 17),
       col = c(oi["red"], oi["blue"], oi["orange"]), cex = 1.3, bty = "n")
plot(0:(nrow(Psd) - 1), err(Psd), type = "l", lwd = 3, col = oi["red"], log = "y",
     ylim = yl, yaxt = "n", xlab = "iteration", ylab = "",
     main = sprintf("steepest descent: %d iterations", nrow(Psd) - 1))
axis(2, at = 10^c(-10, -6, -2, 2), labels = expression(10^-10, 10^-6, 10^-2, 10^2))
dev.off()

# ---- fig05: Nelder-Mead ----------------------------------------------------
svg("fig05-nelder-mead.svg", width = 12, height = 5.2)
par(mfrow = c(1, 2), mar = c(4.6, 4.8, 2.2, 1), cex.lab = 1.5, cex.axis = 1.3,
    las = 1, cex.main = 1.4)
# left: the moves, schematic
S <- rbind(c(0, 0), c(2, 0.4), c(0.6, 1.8))     # worst point = 1st
cen <- colMeans(S[2:3, ])
R <- cen + (cen - S[1, ]); E <- cen + 2 * (cen - S[1, ]); C <- cen - 0.5 * (cen - S[1, ])
plot(NA, xlim = c(-0.3, 4.2), ylim = c(-0.3, 3.6), asp = 1, axes = FALSE, xlab = "", ylab = "",
     main = "moves of the simplex")
polygon(S, border = oi["blue"], lwd = 3)
polygon(rbind(R, S[2, ], S[3, ]), border = oi["orange"], lwd = 2, lty = 2)
points(S, pch = 19, cex = 1.8, col = oi["blue"])
points(rbind(R, E, C), pch = c(17, 15, 18), cex = 2, col = c(oi["orange"], oi["green"], oi["purple"]))
segments(S[1, 1], S[1, 2], E[1], E[2], lty = 3, lwd = 1.5)
text(S[1, 1], S[1, 2] - 0.25, "worst", cex = 1.3)
text(R[1] + 0.35, R[2], "reflect", cex = 1.3, adj = 0, col = oi["orange"])
text(E[1] + 0.2, E[2] - 0.15, "expand", cex = 1.3, adj = 0, col = oi["green"])
text(C[1] + 0.2, C[2] - 0.1, "contract", cex = 1.3, adj = 0, col = oi["purple"])
# right: what optim's Nelder-Mead evaluates on the (EC50, hill) slice
trace <- NULL
f2 <- function(q) { trace <<- rbind(trace, q); sse(c(pstar[1:2], q)) }
o2 <- optim(c(2.6, 0.6), f2)
ec <- seq(0.4, 2.9, length.out = 150); hl <- seq(0.3, 2.0, length.out = 150)
Z <- outer(ec, hl, Vectorize(function(a, b) sse(c(pstar[1:2], a, b))))
contour(ec, hl, Z, levels = c(320, 400, 600, 1000, 2000, 4000, 8000), asp = 1,
        drawlabels = FALSE, col = oi["grey"], lwd = 1.5,
        xlab = expression("EC50 (" * mu * "M)"), ylab = "hill",
        main = sprintf("optim Nelder-Mead: %d evaluations", nrow(trace)))
lines(trace, col = oi["orange"], lwd = 1)
points(trace, pch = 19, cex = 0.8, col = oi["blue"])
points(trace[1, 1], trace[1, 2], pch = 15, cex = 2, col = oi["red"])
points(pstar[3], pstar[4], pch = 4, cex = 2, lwd = 3, col = oi["red"])
cat("NM 2-D slice:", o2$counts, o2$par, "\n")
dev.off()

# ---- fig06: dose in nM, start on a plateau ---------------------------------
open_svg("fig06-units-plateau.svg", 8, 4.8)
dose_nM <- dose * 1000
sse_nM <- function(p) sum((y - f4(p, dose_nM))^2)
ecs <- 10^seq(-1, 6, length.out = 400)
prof <- sapply(ecs, function(e) sse_nM(c(0, 100, e, 1)))
plot(ecs, prof, log = "x", type = "l", lwd = 3, col = oi["blue"], xaxt = "n",
     xlab = "EC50 (nM)", ylab = "SSE", ylim = c(0, max(prof) * 1.05))
axis(1, at = 10^(-1:6), labels = c("0.1", "1", "10", "100", "1000", "1e4", "1e5", "1e6"))
rect(min(dose_nM), -1e6, max(dose_nM), 1e6, col = adjustcolor(oi["sky"], 0.15), border = NA)
points(1, sse_nM(c(0, 100, 1, 1)), pch = 17, cex = 2.2, col = oi["red"])
text(1.4, sse_nM(c(0, 100, 1, 1)) + 3000, "start: EC50 = 1", cex = 1.3, adj = 0, col = oi["red"])
i <- which.min(prof)
points(ecs[i], prof[i], pch = 19, cex = 2, col = oi["orange"])
text(300, 3500, "doses", cex = 1.3, col = oi["blue"])
dev.off()
