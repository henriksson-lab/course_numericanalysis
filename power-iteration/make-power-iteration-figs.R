# Figures for the power-iteration/ unit ("Power iteration and friends").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd power-iteration && Rscript make-power-iteration-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito (colour-blind safe).
# Makes:
#   fig01-stocked-lake.svg   a stocked fish population projected with the Leslie matrix
#   fig02-power-error.svg    power iteration error vs iteration, against |lambda2/lambda1|^k
#   fig03-semelparous.svg    eigenvalues and growth estimates: iteroparous vs semelparous fish
#   fig04-cell-states.svg    cell-state proportions from three sorted starting populations
#   fig05-slow-chain.svg     slow-switching chain: true error, change per step, inverse iteration
#   fig06-rayleigh.svg       eigenvalue error: sum ratio vs Rayleigh quotient, symmetric and not

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")
agecol <- oi[c("blue", "orange", "green", "purple")]

fig <- function(file, w = 9, h = 5.2, expr) {
  svg(file, width = w, height = h)
  par(mar = c(4.8, 5.6, 2.2, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.4)
  expr
  dev.off()
}
logaxis <- function(side, at) {
  axis(side, at = 10^at, labels = parse(text = paste0("10^", at)))
}

## ---- running examples (same numbers as on the slides) ----------------------

L <- matrix(0, 4, 4)
L[1, ] <- c(0, 1.5, 4, 5)                  # offspring per fish of each age
L[cbind(2:4, 1:3)] <- c(0.3, 0.6, 0.5)     # survival to the next age
Ls <- L; Ls[1, ] <- c(0, 0, 0, 12)         # semelparous: spawn once at age 3
evL <- eigen(L)
lam1 <- Re(evL$values[1])
vL <- Re(evL$vectors[, 1]); vL <- vL / sum(vL)
rho <- Mod(evL$values[2]) / lam1

P <- matrix(c(0.80, 0.15, 0.05,
              0.05, 0.90, 0.05,
              0.05, 0.10, 0.85), 3, byrow = TRUE)
Pslow <- diag(3) + 0.01 * (P - diag(3))
piP <- c(0.2, 0.55, 0.25)

C <- matrix(c(1.0, 0.4, 0.3,
              0.4, 1.0, 0.2,
              0.3, 0.2, 1.0), 3)

## ---- fig01: stocked lake ----------------------------------------------------

fig("fig01-stocked-lake.svg", w = 11, h = 4.8, expr = {
  par(mfrow = c(1, 2))
  T <- 30
  N <- matrix(0, T + 1, 4); N[1, ] <- c(0, 1000, 0, 0)
  for (t in 1:T) N[t + 1, ] <- L %*% N[t, ]
  Nl <- N; Nl[Nl == 0] <- NA
  matplot(0:T, Nl, type = "l", lty = 1, lwd = 3, col = agecol, log = "y",
          xlab = "year", ylab = "fish", yaxt = "n", ylim = c(50, 3e5))
  logaxis(2, 2:5)
  legend("topleft", c("age 0", "age 1", "age 2", "age 3"), col = agecol,
         lwd = 3, bty = "n", cex = 1.2)
  title("numbers per age class", font.main = 1)
  prop <- N / rowSums(N)
  matplot(0:T, prop, type = "l", lty = 1, lwd = 3, col = agecol,
          xlab = "year", ylab = "proportion", ylim = c(0, 1))
  abline(h = vL, col = agecol, lty = 2, lwd = 1.5)
  title("age distribution (dashed: eigenvector)", font.main = 1)
})

## ---- fig02: power iteration error -------------------------------------------

fig("fig02-power-error.svg", w = 9, h = 5.2, expr = {
  K <- 75
  x <- c(0, 1, 0, 0); errl <- errx <- numeric(K)
  for (k in 1:K) {
    y <- drop(L %*% x); lam <- sum(y)
    errl[k] <- abs(lam - lam1)
    x <- y / lam
    errx[k] <- max(abs(x - vL))
  }
  errl[errl == 0] <- NA; errx[errx == 0] <- NA
  plot(1:K, errx, log = "y", type = "o", pch = 16, cex = 0.8, lwd = 2,
       col = oi["blue"], ylim = c(1e-16, 10), yaxt = "n",
       xlab = "iteration k", ylab = "error")
  logaxis(2, seq(-16, 0, by = 4))
  lines(1:K, errl, type = "o", pch = 17, cex = 0.8, lwd = 2, col = oi["orange"])
  lines(1:K, 2 * rho^(1:K), lty = 2, lwd = 2.5, col = oi["grey"])
  abline(h = .Machine$double.eps, lty = 3, col = "grey40")
  legend("bottomleft", c("age distribution  max|x - v|",
                         "growth rate  |lambda_k - lambda_1|",
                         expression(2 %*% "|" * lambda[2] / lambda[1] * "|"^k * " = " * 2 %*% 0.654^k)),
         col = oi[c("blue", "orange", "grey")], lwd = 2.5, lty = c(1, 1, 2),
         pch = c(16, 17, NA), bty = "n", cex = 1.2)
})

## ---- fig03: semelparous ------------------------------------------------------

fig("fig03-semelparous.svg", w = 11, h = 4.8, expr = {
  par(mfrow = c(1, 2))
  par(mar = c(4.8, 5.6, 2.2, 1), pty = "s")
  th <- seq(0, 2 * pi, length.out = 200)
  plot(NA, xlim = c(-1.3, 1.3), ylim = c(-1.3, 1.3), asp = 1,
       xlab = "real part", ylab = "imaginary part")
  abline(h = 0, v = 0, col = "grey80")
  lines(lam1 * cos(th), lam1 * sin(th), col = oi["blue"], lty = 2, lwd = 1.5)
  lam_s <- Mod(eigen(Ls)$values[1])
  lines(lam_s * cos(th), lam_s * sin(th), col = oi["red"], lty = 3, lwd = 1.5)
  ev <- evL$values; points(Re(ev), Im(ev), pch = 16, cex = 2, col = oi["blue"])
  ev <- eigen(Ls)$values; points(Re(ev), Im(ev), pch = 17, cex = 2, col = oi["red"])
  title("eigenvalues", font.main = 1)
  par(pty = "m")
  K <- 24
  est <- function(A) {
    x <- c(0, 1, 0, 0); e <- numeric(K)
    for (k in 1:K) { y <- drop(A %*% x); e[k] <- sum(y) / sum(x); x <- y / sum(y) }
    e
  }
  plot(1:K, est(L), type = "o", pch = 16, lwd = 2, col = oi["blue"], log = "y",
       ylim = c(0.2, 80), xlab = "iteration k", ylab = "estimate of growth rate")
  lines(1:K, est(Ls), type = "o", pch = 17, lwd = 2, col = oi["red"])
  abline(h = c(lam1, lam_s), lty = c(2, 3), col = oi[c("blue", "red")])
  legend("top", c("spawns at ages 1-3", "spawns once, at age 3"), ncol = 2,
         pch = c(16, 17), lwd = 2, col = oi[c("blue", "red")], bty = "n", cex = 1.1)
  title("sum(Ax) / sum(x)", font.main = 1)
})

## ---- fig04: cell states -----------------------------------------------------

fig("fig04-cell-states.svg", w = 11, h = 4.2, expr = {
  par(mfrow = c(1, 3), mar = c(4.8, 5, 2.6, 0.8))
  stcol <- oi[c("red", "blue", "green")]
  D <- 30
  for (s in 1:3) {
    p <- matrix(0, D + 1, 3); p[1, s] <- 1
    for (d in 1:D) p[d + 1, ] <- p[d, ] %*% P
    matplot(0:D, p, type = "l", lty = 1, lwd = 3, col = stcol, ylim = c(0, 1),
            xlab = "day", ylab = if (s == 1) "proportion of cells" else "")
    abline(h = piP, col = stcol, lty = 2, lwd = 1.5)
    title(paste("sorted", c("stem-like", "basal", "luminal")[s], "cells"),
          font.main = 1)
    if (s == 1) legend("right", c("S", "B", "L"), col = stcol, lwd = 3,
                       bty = "n", cex = 1.3)
  }
})

## ---- fig05: slow chain --------------------------------------------------------

fig("fig05-slow-chain.svg", w = 9, h = 5.2, expr = {
  K <- 6000
  p <- c(1, 0, 0); err <- ch <- numeric(K)
  for (k in 1:K) {
    pn <- drop(p %*% Pslow); ch[k] <- max(abs(pn - p)); p <- pn
    err[k] <- max(abs(p - piP))
  }
  # inverse iteration, shift 1 + 1e-6
  M <- t(Pslow) - (1 + 1e-6) * diag(3); x <- c(1, 0, 0); ierr <- numeric(3)
  for (k in 1:3) { y <- solve(M, x); x <- y / sum(y); ierr[k] <- max(abs(x - piP)) }
  plot(1:K, err, type = "l", log = "xy", lwd = 3, col = oi["blue"],
       ylim = c(1e-12, 1), xaxt = "n", yaxt = "n",
       xlab = "iteration k", ylab = "")
  axis(1, at = 10^(0:3), labels = c(1, 10, 100, 1000))
  logaxis(2, seq(-12, 0, by = 2))
  lines(1:K, ch, lwd = 3, col = oi["orange"])
  abline(h = 1e-8, lty = 2, col = "grey40")
  text(1.5, 2.2e-8, "tolerance 1e-8", adj = 0, cex = 1.2, col = "grey30")
  k_stop <- which(ch < 1e-8)[1]
  points(k_stop, err[k_stop], pch = 1, cex = 2, lwd = 2)
  points(1:3, ierr, type = "o", pch = 17, cex = 1.6, lwd = 2, col = oi["red"])
  legend("bottomleft", c("power: true error", "power: change per step",
                         "inverse iteration, shift 1 + 1e-6"),
         col = oi[c("blue", "orange", "red")], lwd = 3, pch = c(NA, NA, 17),
         bty = "n", cex = 1.2)
})

## ---- fig06: Rayleigh quotient -------------------------------------------------

fig("fig06-rayleigh.svg", w = 11, h = 4.8, expr = {
  par(mfrow = c(1, 2))
  run <- function(A, x, l1, K) {
    e1 <- e2 <- numeric(K)
    for (k in 1:K) {
      y <- drop(A %*% x)
      e1[k] <- abs(sum(y) / sum(x) - l1)
      e2[k] <- abs(sum(x * y) / sum(x * x) - l1)
      x <- y / sqrt(sum(y^2))
    }
    e1[e1 == 0] <- NA; e2[e2 == 0] <- NA
    cbind(e1, e2)
  }
  panel <- function(E, r, main, squared = TRUE) {
    K <- nrow(E)
    plot(1:K, E[, 1], log = "y", type = "o", pch = 16, cex = 0.8, lwd = 2,
         col = oi["orange"], ylim = c(1e-16, 10), yaxt = "n",
         xlab = "iteration k", ylab = "error in eigenvalue")
    logaxis(2, seq(-16, 0, by = 4))
    lines(1:K, E[, 2], type = "o", pch = 15, cex = 0.8, lwd = 2, col = oi["purple"])
    lines(1:K, r^(1:K), lty = 2, lwd = 2, col = oi["grey"])
    r2 <- r^(2 * (1:K)); r2[r2 < 1e-16] <- NA
    if (squared) lines(1:K, r2, lty = 3, lwd = 2, col = oi["grey"])
    title(main, font.main = 1)
  }
  l1 <- eigen(C)$values[1]
  panel(run(C, c(1, 0, 0), l1, 40), eigen(C)$values[2] / l1,
        "symmetric: trait correlations")
  legend("bottomleft", c("sum(Ax) / sum(x)", "Rayleigh quotient",
                         expression("|" * lambda[2] / lambda[1] * "|"^k),
                         expression("|" * lambda[2] / lambda[1] * "|"^{2 * k})),
         col = oi[c("orange", "purple", "grey", "grey")], lwd = 2,
         lty = c(1, 1, 2, 3), pch = c(16, 15, NA, NA), bty = "n", cex = 1.1)
  panel(run(L, c(0, 1, 0, 0), lam1, 40), rho, "not symmetric: Leslie matrix",
        squared = FALSE)
})
