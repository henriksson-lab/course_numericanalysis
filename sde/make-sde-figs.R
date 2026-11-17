# Figures for the sde/ unit ("Stochastic differential equations").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd sde && Rscript make-sde-figs.R
#
# Base R graphics only; writes fig01-...svg to fig07-...svg (fig05 is PNG:
# many paths). Everything is computed here with the same models and numbers as
# on the slides. Takes about half a minute.

col_a   <- "#0072B2"   # blue
col_b   <- "#D55E00"   # vermillion
col_c   <- "#009E73"   # bluish green
col_d   <- "#CC79A7"   # reddish purple
col_e   <- "#E69F00"   # orange
col_ref <- "grey45"

fig <- function(file, w = 9, h = 5.2, expr, mfrow = c(1, 1)) {
  # PNG for figures with many points (keeps files small), SVG otherwise
  if (grepl("\\.png$", file)) png(file, width = w, height = h, units = "in", res = 150)
  else svg(file, width = w, height = h)
  par(mfrow = mfrow, mar = c(4.8, 5.6, 2.6, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.4)
  expr
  dev.off()
}

## ---- fig01: Brownian motion as the limit of a random walk --------------------

set.seed(3)
nf <- 2^12; hf <- 1 / nf
Wf <- c(0, cumsum(sqrt(hf) * rnorm(nf)))          # one fine path on [0, 1]
tf <- (0:nf) * hf
fig("fig01-random-walk.svg", w = 11, h = 4.8, mfrow = c(1, 2), {
  plot(tf, Wf, type = "n", xlab = "t", ylab = "W(t)", main = "Same path, three step sizes")
  for (R in c(512, 64)) {
    idx <- seq(1, nf + 1, by = R)
    lines(tf[idx], Wf[idx], col = if (R == 512) col_b else col_c, lwd = 2.2, type = "o",
          pch = 16, cex = if (R == 512) 0.9 else 0.4)
  }
  lines(tf, Wf, col = col_a, lwd = 0.8)
  legend("bottomleft", c("h = 1/8", "h = 1/64", "h = 1/4096"), col = c(col_b, col_c, col_a),
         lwd = c(2.2, 2.2, 1.2), bty = "n", cex = 1.2)
  # zoom: first 1/16 of the time, scaled by 4 in space, looks just as rough
  z <- 1:(nf / 16 + 1)
  plot(tf[z] * 16, Wf[z] * 4, type = "l", col = col_a, lwd = 0.8,
       xlab = "16 t", ylab = "4 W(t)", main = "Zoom in: just as rough")
})

## ---- bead tracks and MSD (slides "A bead in a microchannel", "D from the MSD") --

D <- 0.43; dt <- 0.1                               # µm^2/s, s per frame
track <- function(n, u = 0) {
  xy <- apply(matrix(rnorm(2 * n, sd = sqrt(2 * D * dt)), n), 2, cumsum)
  xy[, 1] <- xy[, 1] + u * dt * (1:n)
  xy
}
msd <- function(xy, K) sapply(1:K, function(k)
  mean(rowSums((xy[-(1:k), , drop = FALSE] - xy[1:(nrow(xy) - k), , drop = FALSE])^2)))
Dhat <- function(xy, K) { lag <- (1:K) * dt; unname(coef(lm(msd(xy, K) ~ 0 + lag)) / 4) }

set.seed(1)
tracks <- lapply(1:20, function(i) track(100))
fig("fig02-tracks.svg", w = 11, h = 5.0, mfrow = c(1, 2), {
  cols <- rep(c(col_a, col_b, col_c, col_d, col_e), 2)
  plot(NA, xlim = c(-4, 4), ylim = c(-4, 4), asp = 1, xlab = "x (µm)", ylab = "y (µm)",
       main = "10 beads, 100 frames (10 s)")
  for (i in 1:10) {
    xy <- rbind(0, tracks[[i]])
    lines(xy[, 1], xy[, 2], col = cols[i], lwd = 1.3)
    points(0, 0, pch = 16)
  }
  lag <- (1:50) * dt
  plot(NA, xlim = c(0, 5), ylim = c(0, 18), xlab = "lag t (s)",
       ylab = expression("MSD (" * mu * m^2 * ")"), main = "MSD of 20 tracks")
  for (i in 1:20) lines(lag, msd(tracks[[i]], 50), col = adjustcolor("grey30", 0.5))
  abline(0, 4 * D, col = col_b, lwd = 3)
  legend("topleft", "4 D t", col = col_b, lwd = 3, bty = "n", cex = 1.3)
})

set.seed(2)
cases <- list(c(20, 4), c(20, 10), c(100, 4), c(100, 50), c(1000, 4), c(1000, 500))
est <- lapply(cases, function(cs) replicate(1000, Dhat(track(cs[1]), cs[2])))
for (i in seq_along(cases))
  cat(sprintf("n %4d  lags 1..%-3d  mean %.3f  sd %.3f\n", cases[[i]][1], cases[[i]][2],
              mean(est[[i]]), sd(est[[i]])))
fig("fig03-msd-fit.svg", w = 10, h = 5.0, {
  par(mar = c(6.2, 5.6, 2.6, 1))
  boxplot(est, names = FALSE, xaxt = "n",
          col = rep(c(col_c, col_d), 3), outline = FALSE, ylim = c(0, 1.2),
          ylab = expression(hat(D) ~ "(" * mu * m^2 / s * ")"), xlab = "",
          main = "1000 tracks each: D from a line through MSD(lags 1..K)")
  axis(1, at = 1:6, tick = FALSE, line = 0.3, cex.axis = 1.25,
       labels = sapply(cases, function(cs) paste0(cs[1], " frames")))
  axis(1, at = 1:6, tick = FALSE, line = 1.9, cex.axis = 1.25,
       labels = sapply(cases, function(cs) paste0("K = ", cs[2])))
  abline(h = D, col = col_b, lwd = 2, lty = 2)
})

## ---- fig04: strong and weak error of EM and Milstein on GBM ------------------

r <- 0.8; s <- 0.5; N0 <- 1; T <- 1
set.seed(1); M <- 10000; nf <- 2^10; hf <- T / nf
dWf <- matrix(sqrt(hf) * rnorm(nf * M), nf)
NT <- N0 * exp((r - s^2 / 2) * T + s * colSums(dWf))
step <- function(dW, h, milstein = FALSE) {
  N <- rep(N0, ncol(dW))
  for (i in 1:nrow(dW)) N <- N + r * N * h + s * N * dW[i, ] +
    if (milstein) 0.5 * s^2 * N * (dW[i, ]^2 - h) else 0
  N
}
Rs <- 2^(0:8); hs <- Rs * hf
err <- t(sapply(Rs, function(R) {
  dW <- rowsum(dWf, rep(1:(nf / R), each = R))
  Ne <- step(dW, R * hf); Nm <- step(dW, R * hf, TRUE)
  c(strong = mean(abs(Ne - NT)), weak = abs(mean(Ne - NT)), mil = mean(abs(Nm - NT)))
}))
print(cbind(h = hs, err))
cat("slopes:", coef(lm(log10(err[, "strong"]) ~ log10(hs)))[2],
    coef(lm(log10(err[, "weak"]) ~ log10(hs)))[2], coef(lm(log10(err[, "mil"]) ~ log10(hs)))[2], "\n")
fig("fig04-convergence.svg", w = 8.5, h = 5.4, {
  par(mar = c(4.8, 5.6, 1, 1))
  plot(hs, err[, "strong"], log = "xy", type = "o", pch = 16, lwd = 2.5, col = col_a,
       ylim = c(3e-4, 0.5), xlab = "step h (hours)", ylab = "error at T = 1 h", xaxt = "n", yaxt = "n")
  axis(1, at = 2^-(10:2), labels = paste0("1/", 2^(10:2)))
  axis(2, at = 10^(-3:0), labels = c("0.001", "0.01", "0.1", "1"))
  lines(hs, err[, "weak"], type = "o", pch = 17, lwd = 2.5, col = col_b)
  lines(hs, err[, "mil"], type = "o", pch = 15, lwd = 2.5, col = col_c)
  lines(hs, 0.3 * sqrt(hs), lty = 2, col = col_ref, lwd = 1.5)
  lines(hs, 0.8 * hs, lty = 3, col = col_ref, lwd = 1.5)
  abline(h = 2 * 1.186 / sqrt(M), col = col_d, lty = 4, lwd = 2)
  legend("bottomright", c("EM strong (slope 0.5)", "EM weak (slope 1)", "Milstein strong (slope 1)",
                          "2 SE of a plain mean, M = 10 000"),
         col = c(col_a, col_b, col_c, col_d), lwd = 2.5, pch = c(16, 17, 15, NA), lty = c(1, 1, 1, 4),
         bty = "n", cex = 1.1)
})
rm(dWf)

## ---- fig05: GBM, mean vs median (Ito) -----------------------------------------

r <- 0.1; s <- 0.6; N0 <- 1000; h <- 0.05; Tg <- 50
set.seed(1); P <- 10000; nst <- Tg / h
N <- matrix(0, nst + 1, 100); Nall <- rep(N0, P); N[1, ] <- N0
med <- numeric(nst + 1); med[1] <- N0; low <- numeric(nst + 1)
for (i in 1:nst) {
  Nall <- Nall + r * Nall * h + s * Nall * sqrt(h) * rnorm(P)
  N[i + 1, ] <- Nall[1:100]; med[i + 1] <- median(Nall); low[i + 1] <- mean(Nall < 10)
}
tg <- (0:nst) * h
cat("GBM: median(50) =", med[nst + 1], " P(N(50) < 10) =", low[nst + 1], "\n")
fig("fig05-gbm.png", w = 9, h = 5.2, {
  par(mar = c(4.8, 5.6, 1, 1))
  matplot(tg, pmax(N, 1e-3), type = "l", lty = 1, log = "y", col = adjustcolor("grey30", 0.25),
          xlab = "time (years)", ylab = "population N", ylim = c(1e-1, 1e6), yaxt = "n")
  axis(2, at = 10^(-1:6), labels = c("0.1", "1", "10", "100", "1000", expression(10^4),
                                     expression(10^5), expression(10^6)))
  lines(tg, N0 * exp(r * tg), col = col_b, lwd = 3.5)
  lines(tg, N0 * exp((r - s^2 / 2) * tg), col = col_a, lwd = 3.5)
  lines(tg, med, col = col_c, lwd = 2, lty = 2)
  abline(h = 10, lty = 3)
  legend("topleft", c(expression("mean " * e^{r * t}), expression("median " * e^{(r - sigma^2/2) * t}),
                      "median of 10 000 EM paths", "100 single paths"),
         col = c(col_b, col_a, col_c, "grey50"), lwd = c(3.5, 3.5, 2, 1), lty = c(1, 1, 2, 1),
         bty = "n", cex = 1.15)
})

## ---- fig06: mean-square stability, EM vs drift-implicit EM ---------------------

k <- 0.8; mu <- 1; C0 <- 10; Th <- 48
fig("fig06-ms-stability.svg", w = 9, h = 5.2, {
  par(mar = c(4.8, 5.6, 1, 1))
  # E[C_n^2] of each method follows from its one-step factor (exercise B5)
  plot(NA, xlim = c(0, Th), ylim = c(1e-12, 1e12), log = "y", xlab = "time (hours)",
       ylab = expression("mean square  E[" * C[n]^2 * "]"), yaxt = "n")
  axis(2, at = 10^seq(-12, 12, 6), labels = c(expression(10^-12), expression(10^-6), "1",
                                              expression(10^6), expression(10^12)))
  meth <- function(h, implicit, col, pch) {
    fac <- if (implicit) (1 + mu^2 * h) / (1 + k * h)^2 else (1 - k * h)^2 + mu^2 * h
    tt <- (0:(Th / h)) * h
    lines(tt, C0^2 * fac^(tt / h), col = col, lwd = 2.5, type = "o", pch = pch, cex = 0.8)
  }
  tt <- seq(0, Th, 0.5)
  lines(tt, C0^2 * exp((-2 * k + mu^2) * tt), lwd = 3)
  meth(0.5, FALSE, col_a, 16); meth(1, FALSE, col_e, 18); meth(2, FALSE, col_b, 17)
  meth(2, TRUE, col_c, 15)
  legend("topleft", c("exact SDE", "EM, h = 0.5 h", "EM, h = 1 h", "EM, h = 2 h",
                      "drift-implicit EM, h = 2 h"),
         col = c("black", col_a, col_e, col_b, col_c), lwd = c(3, 2.5, 2.5, 2.5, 2.5),
         pch = c(NA, 16, 18, 17, 15), bty = "n", cex = 1.15)
})

## ---- fig07: Fokker-Planck PDE vs histogram of EM paths --------------------------

theta <- 1; D <- 430; sigma <- sqrt(2 * D)
dx <- 2; x <- seq(-100, 160, by = dx); xf <- x[-1] - dx / 2; u <- -theta * xf
dtp <- 0.002; p <- dnorm(x, 60, 5)
set.seed(1); X <- rnorm(10000, 60, 5); hE <- 0.01
snap <- c(0.2, 1, 3); Pk <- list(); Xk <- list(); tnow <- 0
for (j in seq_along(snap)) {
  while (tnow < snap[j] - 1e-9) {
    for (q in 1:round(hE / dtp)) {
      F <- u * (head(p, -1) + p[-1]) / 2 - D * diff(p) / dx
      p <- p - dtp / dx * (c(F, 0) - c(0, F))
    }
    X <- X - theta * X * hE + sigma * sqrt(hE) * rnorm(10000)
    tnow <- tnow + hE
  }
  Pk[[j]] <- p; Xk[[j]] <- X
}
fig("fig07-fokker-planck.svg", w = 11, h = 4.8, mfrow = c(1, 3), {
  par(cex = 1)
  for (j in seq_along(snap)) {
    par(mar = c(4.8, if (j == 1) 5.6 else 2, 2.6, 0.6))
    hh <- hist(Xk[[j]], breaks = seq(-130, 190, 8), plot = FALSE)
    plot(hh, freq = FALSE, col = "grey85", border = "grey60", xlim = c(-80, 140), ylim = c(0, 0.04),
         main = paste0("t = ", snap[j], " ms"), xlab = "bead position x (nm)",
         ylab = if (j == 1) "density" else "", yaxt = if (j == 1) "s" else "n")
    lines(x, Pk[[j]], col = col_a, lwd = 3)
    ex_m <- 60 * exp(-theta * snap[j])
    ex_v <- D / theta * (1 - exp(-2 * theta * snap[j])) + 25 * exp(-2 * theta * snap[j])
    lines(x, dnorm(x, ex_m, sqrt(ex_v)), col = col_b, lwd = 2, lty = 2)
    if (j == 1) legend("topleft", c("10 000 EM paths", "Fokker-Planck PDE", "exact"),
                       col = c("grey60", col_a, col_b), lwd = c(8, 3, 2), lty = c(1, 1, 2),
                       bty = "n", cex = 1.0)
  }
})
