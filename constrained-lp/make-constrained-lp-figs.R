# Figures for the constrained-lp/ unit ("Constraints and linear programming").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd constrained-lp && Rscript make-constrained-lp-figs.R
#
# Base R graphics only (grDevices::svg) plus lpSolve. Colours are Okabe-Ito
# (colour-blind safe). Makes:
#   fig01-decay-bounds.svg   mRNA decay with background: free fit vs bounded fit
#   fig02-lagrange.svg       Lagrange picture: contours, constraint, gradients
#   fig03-lp-feasible.svg    the ATP LP: feasible polygon, objective lines
#   fig04-simplex-path.svg   simplex walk (0,0) -> (3,0) -> (2,4)
#   fig05-lp-failures.svg    degenerate / unbounded / infeasible
#   fig06-network.svg        the toy metabolic network (9 reactions, 6 metabolites)
#   fig07-fva.svg            flux variability ranges at maximal growth
#   fig08-oxygen-scan.svg    growth, acetate and O2 shadow price vs O2 uptake limit

library(lpSolve)

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- fig01: decay with a background that must be >= 0 ---------------------
t <- 0:8
y <- c(1.020, 0.577, 0.378, 0.243, 0.125, 0.062, 0.020, 0.010, -0.012)
sse <- function(p) sum((y - (p[1] * exp(-p[2] * t) + p[3]))^2)
free <- optim(c(1, 0.3, 0.1), sse, method = "BFGS")$par
bnd  <- optim(c(1, 0.3, 0.1), sse, method = "L-BFGS-B", lower = c(0, 0, 0))$par
tt <- seq(0, 8, length.out = 200)

svg("fig01-decay-bounds.svg", width = 10, height = 5.2)
par(mar = c(4.8, 5.2, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(t, y, pch = 19, cex = 1.6, ylim = c(-0.08, 1.05),
     xlab = "time after transcription block (h)", ylab = "mRNA signal")
abline(h = 0, col = oi["grey"])
lines(tt, free[1] * exp(-free[2] * tt) + free[3], col = oi["red"], lwd = 3, lty = 2)
lines(tt, bnd[1] * exp(-bnd[2] * tt) + bnd[3], col = oi["blue"], lwd = 3)
legend("topright", bty = "n", cex = 1.35, lwd = 3, lty = c(2, 1),
       col = oi[c("red", "blue")],
       legend = c(sprintf("free fit: background c = %.3f", free[3]),
                  sprintf("bound c >= 0: c = %.3f", bnd[3])))
dev.off()

## ---- fig02: Lagrange multipliers in one picture ---------------------------
## maximise f = 30 log p1 + 10 log p2 subject to p1 + p2 = 1
f <- function(p1, p2) 30 * log(p1) + 10 * log(p2)
g1 <- seq(0.05, 1.1, length.out = 300)
z <- outer(g1, g1, f)
svg("fig02-lagrange.svg", width = 7, height = 6.4)
par(mar = c(4.8, 5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
lev <- f(0.75, 0.25) + c(-12, -6, -3, 0, 3, 6, 10)
contour(g1, g1, z, levels = lev, drawlabels = FALSE, col = oi["grey"], lwd = 1.5,
        xlab = expression(p[1]), ylab = expression(p[2]), asp = 1,
        xlim = c(0.05, 1.1), ylim = c(0.05, 1.1))
contour(g1, g1, z, levels = f(0.75, 0.25), drawlabels = FALSE, add = TRUE,
        col = oi["sky"], lwd = 3)
abline(1, -1, col = oi["blue"], lwd = 4)
text(0.18, 0.93, expression(p[1] + p[2] == 1), col = oi["blue"], cex = 1.4, adj = 0)
# optimum: grad f = (40, 40) = lambda * (1, 1), lambda = 40
s <- 0.006
arrows(0.75, 0.25, 0.75 + s * 40, 0.25 + s * 40, col = oi["red"], lwd = 4, length = 0.12)
arrows(0.75, 0.25, 0.75 + 0.12, 0.25 + 0.12, col = oi["blue"], lwd = 4, length = 0.12)
points(0.75, 0.25, pch = 19, cex = 2)
text(0.75 + s * 40, 0.25 + s * 40 + 0.05, expression(nabla * f == 40 * nabla * g),
     col = oi["red"], cex = 1.3)
# a non-optimal point on the line: grad f = (30/0.4, 10/0.6)
p <- c(0.4, 0.6); gf <- c(30 / p[1], 10 / p[2])
arrows(p[1], p[2], p[1] + 0.003 * gf[1], p[2] + 0.003 * gf[2], col = oi["orange"],
       lwd = 4, length = 0.12)
points(p[1], p[2], pch = 19, cex = 1.6, col = oi["orange"])
text(0.02, 0.38, "not parallel:\nslide along\nthe line", col = oi["orange"],
     cex = 1.25, adj = 0)
text(0.93, 0.31, expression(nabla * g), col = oi["blue"], cex = 1.3)
dev.off()

## ---- the ATP LP: max 6 x1 + 3 x2, x1 + x2 <= 6, 4 x1 + x2 <= 12 ----------
lp_base <- function(xlim = c(0, 7), ylim = c(0, 13)) {
  plot(NA, xlim = xlim, ylim = ylim, xaxs = "i", yaxs = "i",
       xlab = expression("respiration flux " * x[1]),
       ylab = expression("fermentation flux " * x[2]))
}
svg("fig03-lp-feasible.svg", width = 7.5, height = 6.4)
par(mar = c(4.8, 5, 1, 1), cex.lab = 1.4, cex.axis = 1.3)
lp_base()
polygon(c(0, 3, 2, 0), c(0, 0, 4, 6), col = adjustcolor(oi["sky"], 0.35), border = NA)
abline(6, -1, col = oi["blue"], lwd = 3)      # carbon
abline(12, -4, col = oi["green"], lwd = 3)    # proteome
for (zz in c(9, 18, 24, 30)) abline(zz / 3, -2, col = oi["orange"], lty = 2, lwd = 2)
points(c(0, 3, 2, 0), c(0, 0, 4, 6), pch = 19, cex = 1.5)
points(2, 4, pch = 19, cex = 2.6, col = oi["red"])
legend(2.1, 5.3, "optimum (2, 4)\nATP = 24", bty = "o", bg = "white", box.col = NA,
       cex = 1.3, text.col = oi["red"])
text(4.6, 1.0, "carbon\nx1 + x2 <= 6", col = oi["blue"], cex = 1.2)
text(2.4, 9.5, "proteome\n4 x1 + x2 <= 12", col = oi["green"], cex = 1.2)
text(5.3, 11.3, "ATP = 6 x1 + 3 x2\n= 30 (too high)", col = oi["orange"], cex = 1.2)
dev.off()

svg("fig04-simplex-path.svg", width = 7.5, height = 6.4)
par(mar = c(4.8, 5, 1, 1), cex.lab = 1.4, cex.axis = 1.3)
lp_base(ylim = c(0, 7))
polygon(c(0, 3, 2, 0), c(0, 0, 4, 6), col = adjustcolor(oi["sky"], 0.35),
        border = oi["blue"], lwd = 2)
arrows(0.05, 0.12, 2.9, 0.12, col = oi["red"], lwd = 5, length = 0.15)
arrows(3, 0.1, 2.08, 3.75, col = oi["red"], lwd = 5, length = 0.15)
points(c(0, 3, 2, 0), c(0, 0, 4, 6), pch = 19, cex = 1.6)
text(0.1, 0.55, "start (0,0)\nATP = 0", adj = 0, cex = 1.2)
text(3.1, 0.5, "(3,0)\nATP = 18", adj = 0, cex = 1.2)
text(2.15, 4.3, "(2,4)  ATP = 24\nno better neighbour", adj = 0, cex = 1.2,
     col = oi["red"])
text(0.1, 6.3, "(0,6)  ATP = 18", adj = 0, cex = 1.2)
dev.off()

## ---- fig05: three ways an LP goes wrong -----------------------------------
svg("fig05-lp-failures.svg", width = 13, height = 4.6)
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1), cex.lab = 1.4, cex.axis = 1.2,
    cex.main = 1.6)
# (a) tie: objective parallel to the proteome edge: max 8 x1 + 2 x2
lp_base(ylim = c(0, 7)); title("objective parallel to an edge")
polygon(c(0, 3, 2, 0), c(0, 0, 4, 6), col = adjustcolor(oi["sky"], 0.35), border = NA)
segments(3, 0, 2, 4, col = oi["red"], lwd = 6)
for (zz in c(12, 24)) abline(zz / 2, -4, col = oi["orange"], lty = 2, lwd = 2)
text(3.1, 2.6, "every point on the\nred edge is optimal", adj = 0, cex = 1.2, col = oi["red"])
# (b) unbounded: no carbon limit and fermentation costs no protein
lp_base(ylim = c(0, 13)); title("unbounded")
polygon(c(0, 3, 3, 0), c(0, 0, 13, 13), col = adjustcolor(oi["sky"], 0.35), border = NA)
abline(v = 3, col = oi["green"], lwd = 3)
arrows(1.5, 6, 1.5, 12.3, col = oi["red"], lwd = 5, length = 0.15)
text(3.2, 6, "ATP grows\nwithout limit", adj = 0, cex = 1.2, col = oi["red"])
# (c) infeasible: demand x1 + x2 >= 8 but carbon <= 6
lp_base(xlim = c(0, 9), ylim = c(0, 13)); title("infeasible")
polygon(c(0, 6, 0), c(0, 0, 6), col = adjustcolor(oi["sky"], 0.35), border = NA)
polygon(c(8, 9, 9, 0, 0), c(0, 0, 13, 13, 8), col = adjustcolor(oi["orange"], 0.3),
        border = NA)
abline(6, -1, col = oi["blue"], lwd = 3); abline(8, -1, col = oi["orange"], lwd = 3)
text(0.3, 2, "carbon <= 6", adj = 0, cex = 1.2, col = oi["blue"])
text(4.2, 9.5, "demand >= 8", adj = 0, cex = 1.2, col = oi["orange"])
dev.off()

## ---- the toy network -------------------------------------------------------
rxn <- c("EX_glc", "EX_o2", "GLY", "BYP1", "BYP2", "RES", "FER", "EX_ac", "BIO")
S <- rbind(
  G   = c(1, 0, -1, -1,  0,  0,  0,  0,  0),
  O2  = c(0, 1,  0,  0,  0, -2,  0,  0,  0),
  X   = c(0, 0,  0,  1, -1,  0,  0,  0,  0),
  P   = c(0, 0,  2,  0,  2, -1, -1,  0, -1),
  ATP = c(0, 0,  2,  0,  2, 11,  1,  0, -5),
  Ac  = c(0, 0,  0,  0,  0,  0,  1, -1,  0))
colnames(S) <- rxn
n <- ncol(S)
fba <- function(ub = c(10, 15), fix_bio = NULL, obj = c(rep(0, n - 1), 1),
                dir = "max") {
  A <- rbind(S, diag(n)[1:2, ]); dirs <- c(rep("=", 6), "<=", "<=")
  rhs <- c(rep(0, 6), ub)
  if (!is.null(fix_bio)) { A <- rbind(A, diag(n)[n, ]); dirs <- c(dirs, "=")
                           rhs <- c(rhs, fix_bio) }
  lp(dir, obj, A, dirs, rhs, compute.sens = TRUE)
}

svg("fig06-network.svg", width = 12, height = 6)
par(mar = c(0, 0, 0, 0))
plot(NA, xlim = c(-0.3, 11.6), ylim = c(-0.2, 7.6), axes = FALSE, xlab = "", ylab = "")
met <- list(G = c(2, 4.5), X = c(3.75, 2.3), P = c(5.5, 4.5), O2 = c(7.6, 6.6),
            Ac = c(7.6, 0.6))
rx <- list(EX_glc = c(0.4, 4.5), GLY = c(3.75, 5.6), BYP1 = c(2.6, 3.0),
           BYP2 = c(4.9, 3.0), RES = c(7.6, 4.5), FER = c(6.3, 2.2),
           EX_ac = c(9.6, 0.6), EX_o2 = c(7.6, 7.5), BIO = c(10.2, 4.5))
edge <- function(a, b, lab = "", col = "black", lty = 1, lwd = 2.5, sh = 0.42) {
  d <- b - a; L <- sqrt(sum(d^2)); u <- d / L
  p0 <- a + u * sh; p1 <- b - u * sh
  arrows(p0[1], p0[2], p1[1], p1[2], length = 0.12, lwd = lwd, col = col, lty = lty)
  if (nzchar(lab)) text(mean(c(p0[1], p1[1])) + 0.12, mean(c(p0[2], p1[2])) + 0.18,
                        lab, cex = 1.2, col = col)
}
M <- function(m) met[[m]]; R <- function(r) rx[[r]]
edge(R("EX_glc"), M("G"))
edge(M("G"), R("GLY")); edge(R("GLY"), M("P"), "2")
edge(M("G"), R("BYP1")); edge(R("BYP1"), M("X")); edge(M("X"), R("BYP2"))
edge(R("BYP2"), M("P"), "2")
edge(M("P"), R("RES")); edge(R("EX_o2"), M("O2"), sh = 0.35); edge(M("O2"), R("RES"), "2")
edge(M("P"), R("FER")); edge(R("FER"), M("Ac")); edge(M("Ac"), R("EX_ac"))
edge(M("P"), R("BIO"))
ac <- oi["purple"]
atp <- list(GLY = c("+2 ATP", 0, 0.45), BYP2 = c("+2 ATP", 0.2, -0.45),
            RES = c("+11 ATP", 0.95, 0.35), FER = c("+1 ATP", 0.95, 0),
            BIO = c("-5 ATP", 0, -0.45))
for (r in names(atp)) text(rx[[r]][1] + as.numeric(atp[[r]][2]),
                           rx[[r]][2] + as.numeric(atp[[r]][3]), atp[[r]][1],
                           col = ac, cex = 1.2, font = 2)
for (m in names(met)) {
  symbols(met[[m]][1], met[[m]][2], circles = 0.36, inches = FALSE, add = TRUE,
          bg = adjustcolor(oi["sky"], 0.5),
          fg = "black")
  text(met[[m]][1], met[[m]][2], m, cex = 1.3, font = 2)
}
for (r in names(rx)) {
  w <- 0.14 * nchar(r) + 0.15
  col <- if (r %in% c("EX_glc", "EX_o2")) oi["green"] else
         if (r == "BIO") oi["red"] else "white"
  rect(rx[[r]][1] - w / 2, rx[[r]][2] - 0.22, rx[[r]][1] + w / 2, rx[[r]][2] + 0.22,
       col = adjustcolor(col, 0.5), border = "black")
  text(rx[[r]][1], rx[[r]][2], r, cex = 1.1)
}
text(0.4, 5.1, "glucose\n<= 10", cex = 1.1, col = "black")
text(8.75, 7.5, "oxygen <= 15", cex = 1.1, adj = 0)
text(10.2, 5.1, "growth", cex = 1.2, col = oi["red"], font = 2)
text(9.6, 0.05, "acetate out", cex = 1.1)
text(0, 1.2, "purple: ATP made (+) or used (-)", cex = 1.2, col = ac, adj = 0)
dev.off()

## ---- fig07: flux variability at maximal growth ---------------------------
opt <- fba()$objval
rng <- t(sapply(seq_len(n), function(j) {
  o <- numeric(n); o[j] <- 1
  c(fba(fix_bio = opt, obj = o, dir = "min")$objval,
    fba(fix_bio = opt, obj = o, dir = "max")$objval)
}))
svg("fig07-fva.svg", width = 10, height = 5.4)
par(mar = c(4.8, 7.5, 1, 1), cex.lab = 1.4, cex.axis = 1.25)
plot(NA, xlim = c(0, 16), ylim = c(0.5, n + 0.5), yaxt = "n",
     xlab = "flux at maximal growth (min to max)", ylab = "")
axis(2, at = n:1, labels = rxn, las = 1)
fixed <- abs(rng[, 2] - rng[, 1]) < 1e-9
for (j in seq_len(n)) {
  yy <- n + 1 - j
  if (fixed[j]) points(rng[j, 1], yy, pch = 19, cex = 1.8, col = oi["blue"])
  else segments(rng[j, 1], yy, rng[j, 2], yy, lwd = 10, col = oi["orange"], lend = 1)
}
legend("topright", bty = "n", cex = 1.3, pch = c(19, NA), lwd = c(NA, 10),
       col = oi[c("blue", "orange")], legend = c("unique value", "range"))
dev.off()

## ---- fig08: scan the oxygen limit ------------------------------------------
o2 <- seq(0, 15, by = 0.25)
sc <- t(sapply(o2, function(u) { r <- fba(ub = c(10, u))
  c(r$objval, r$solution[8], r$duals[8]) }))
svg("fig08-oxygen-scan.svg", width = 10, height = 5.2)
par(mar = c(4.8, 5.2, 1, 1), cex.lab = 1.45, cex.axis = 1.3)
plot(o2, sc[, 1], type = "l", lwd = 4, col = oi["red"], ylim = c(0, 16),
     xlab = "oxygen uptake limit", ylab = "flux")
lines(o2, sc[, 2], lwd = 4, col = oi["purple"], lty = 2)
abline(v = 10, col = oi["grey"], lty = 3, lwd = 2)
text(5, 13.7, "O2 limits growth\nshadow price 0.833", cex = 1.25)
text(12.6, 11.5, "O2 in excess\nshadow price 0", cex = 1.25)
legend("bottomright", bty = "n", cex = 1.3, lwd = 4, lty = c(1, 2),
       col = oi[c("red", "purple")], legend = c("growth (BIO)", "acetate out (EX_ac)"))
dev.off()
