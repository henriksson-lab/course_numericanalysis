# Figures for the gillespie/ unit ("Stochastic simulation").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd gillespie && Rscript make-gillespie-figs.R
#
# Base R graphics only; writes fig01-...svg to fig07-...svg. Everything is
# computed here, with the same code, seeds and numbers as on the slides.
# Takes about a minute (hand-written SSA loops in R).

col_a   <- "#0072B2"   # blue
col_b   <- "#D55E00"   # vermillion
col_c   <- "#009E73"   # bluish green
col_ref <- "grey45"
col_pur <- "#CC79A7"   # reddish purple
col_run <- adjustcolor("grey30", 0.55)

fig <- function(file, w = 9, h = 5.2, expr, mfrow = c(1, 1)) {
  # PNG for figures with many points (keeps files small), SVG otherwise
  if (grepl("\\.png$", file)) png(file, width = w, height = h, units = "in", res = 150)
  else svg(file, width = w, height = h)
  par(mfrow = mfrow, mar = c(4.8, 5.6, 2.6, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.4)
  expr
  dev.off()
}

## ---- the SSA as on the slides ------------------------------------------------

ssa <- function(x, S, propensity, times) {   # S: species x reactions
  out <- matrix(0, length(times), length(x)); t <- 0; i <- 1
  repeat {
    a <- propensity(x); a0 <- sum(a)
    t <- t - log(runif(1)) / a0               # Inf if a0 = 0: nothing more happens
    while (i <= length(times) && times[i] < t) { out[i, ] <- x; i <- i + 1 }
    if (i > length(times)) return(out)
    j <- findInterval(runif(1) * a0, cumsum(a)) + 1   # reaction j w.p. a_j / a0
    x <- x + S[, j]
  }
}
# the same, but recording every event (for plotting trajectories)
ssa_path <- function(x, S, propensity, t_end) {
  n_max <- 1e6; tt <- numeric(n_max); X <- matrix(0, n_max, length(x))
  t <- 0; n <- 1; X[1, ] <- x
  repeat {
    a <- propensity(x); a0 <- sum(a)
    t <- t - log(runif(1)) / a0
    if (t > t_end) break
    j <- findInterval(runif(1) * a0, cumsum(a)) + 1
    x <- x + S[, j]
    n <- n + 1; tt[n] <- t; X[n, ] <- x
  }
  n <- n + 1; tt[n] <- t_end; X[n, ] <- x      # hold the last state to t_end
  list(t = tt[1:n], x = X[1:n, , drop = FALSE])
}
fano <- function(v) var(v) / mean(v)

S1 <- matrix(c(1, -1), nrow = 1)
bd <- function(k) function(x) c(k, 1 * x)
S3 <- rbind(g = c(1, -1, 0,  0, 0,  0),
            m = c(0,  0, 1, -1, 0,  0),
            p = c(0,  0, 0,  0, 1, -1))
telegraph <- function(x) c(0.5 * (1 - x[1]), 4.5 * x[1], 100 * x[1],
                           1 * x[2], 2 * x[2], 0.2 * x[3])
constit   <- function(x) c(0, 0, 10, 1 * x[2], 2 * x[2], 0.2 * x[3])

## ---- 1000 cells of each gene (slide "Telegraph in R") ------------------------

set.seed(1)
Xc <- t(replicate(1000, ssa(c(1, 0, 0), S3, constit, 30)[1, ]))
set.seed(1)
Xb <- t(replicate(1000, ssa(c(0, 0, 0), S3, telegraph, 30)[1, ]))

## ---- fig01: two genes, same mean, different spread ---------------------------

fig("fig01-two-genes.svg", w = 11, h = 5.0, mfrow = c(1, 2), {
  br <- seq(-0.5, 80.5, 1)
  for (k in 1:2) {
    v <- if (k == 1) Xc[, 2] else Xb[, 2]
    h <- hist(pmin(v, 80), breaks = br, plot = FALSE)
    plot(h$mids, h$counts / 1000, type = "h", lwd = 4, lend = 1,
         col = if (k == 1) col_a else col_b, xlim = c(0, 60), ylim = c(0, 0.2),
         xlab = "mRNA molecules per cell", ylab = "fraction of cells",
         main = if (k == 1) "gene A: constitutive" else "gene B: bursty")
    abline(v = 10, lty = 2, col = col_ref, lwd = 2)
    legend("topright", bty = "n", cex = 1.3,
           legend = c(sprintf("mean %.1f", mean(v)), sprintf("Fano %.1f", fano(v)),
                      sprintf("%.0f%% cells with 0", 100 * mean(v == 0))))
  }
})

## ---- fig02: small and large copy numbers -------------------------------------

set.seed(2)
fig("fig02-small-numbers.png", w = 11, h = 5.0, mfrow = c(1, 2), {
  for (k in c(10, 1000)) {
    plot(NA, xlim = c(0, 10), ylim = c(0, 1.8), xlab = "time (h)",
         ylab = "molecules / steady state",
         main = sprintf("k = %g per hour: mean %g molecules", k, k))
    for (r in 1:3) {
      p <- ssa_path(0, S1, bd(k), 10)
      lines(p$t, p$x[, 1] / k, type = "s", col = c(col_b, col_c, col_pur)[r], lwd = 1.6)
    }
    curve(1 - exp(-x), add = TRUE, lwd = 4, col = col_a)
    legend("topleft", bty = "n", cex = 1.2, lwd = c(4, 1.6), col = c(col_a, col_b),
           legend = c("ODE", "3 simulated cells"))
  }
})

## ---- fig03: one SSA step (m = 10: a = (10, 10)) ------------------------------

fig("fig03-ssa-step.svg", w = 11, h = 4.6, mfrow = c(1, 2), {
  a0 <- 20
  curve(1 - exp(-a0 * x), 0, 0.2, lwd = 4, col = col_a, ylim = c(0, 1),
        xlab = expression("waiting time " * tau * " (h)"),
        ylab = expression(F(tau) == 1 - e^{-a[0] * tau}),
        main = expression("1. when?  " * tau * " = " * -log(u[1]) / a[0]))
  u1 <- 0.3; tau <- -log(u1) / a0
  # inverse transform of 1 - u: plot F and read off at 1 - u1
  segments(0, 1 - u1, tau, 1 - u1, col = col_b, lwd = 2.5)
  arrows(tau, 1 - u1, tau, 0, col = col_b, lwd = 2.5, length = 0.12)
  text(tau, 0.08, sprintf("  %.4f h", tau), adj = 0, cex = 1.3, col = col_b)
  text(0.003, 1 - u1 + 0.07, expression(1 - u[1] == 0.7), adj = 0, cex = 1.3, col = col_b)
  par(mar = c(4.8, 2, 2.6, 1))
  plot(NA, xlim = c(-0.5, 21.5), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "",
       main = expression("2. which?  " * u[2] * a[0] * " falls in reaction j's piece"))
  rect(0, 0.35, 10, 0.65, col = adjustcolor(col_a, 0.35), border = col_a, lwd = 2)
  rect(10, 0.35, 20, 0.65, col = adjustcolor(col_c, 0.35), border = col_c, lwd = 2)
  text(5, 0.5, expression("make: " * a[1] == 10), cex = 1.3)
  text(15, 0.5, expression("degrade: " * a[2] == 10), cex = 1.3)
  axis(1, at = c(0, 10, 20), labels = c("0", expression(a[1]), expression(a[0] == 20)))
  arrows(12.4, 0.95, 12.4, 0.67, lwd = 2.5, col = col_b, length = 0.12)
  text(12.4, 0.98, expression(u[2] * a[0] == "0.62 x 20 = 12.4"), cex = 1.3, col = col_b)
  text(10, 0.15, "degrade: m goes from 10 to 9", cex = 1.3)
})

## ---- fig04: many runs vs ODE, distribution vs Poisson ------------------------

set.seed(1)
tm <- c(0.5, 1, 2, 10)
M <- replicate(2000, ssa(0, S1, bd(10), tm)[, 1])        # as on the slide
set.seed(3)
tg <- seq(0, 10, by = 0.25)
Mg <- replicate(2000, ssa(0, S1, bd(10), tg)[, 1])
fig("fig04-bd-runs.svg", w = 11, h = 5.0, mfrow = c(1, 2), {
  plot(NA, xlim = c(0, 10), ylim = c(0, 22), xlab = "time (h)", ylab = "mRNA molecules",
       main = "5 cells, mean of 2000, ODE")
  for (r in 1:5) { p <- ssa_path(0, S1, bd(10), 10); lines(p$t, p$x[, 1], type = "s", col = col_run) }
  curve(10 * (1 - exp(-x)), add = TRUE, lwd = 4, col = col_a)
  points(tg, rowMeans(Mg), pch = 19, col = col_b, cex = 0.9)
  legend("bottomright", bty = "n", cex = 1.2, lwd = c(4, NA, 1), pch = c(NA, 19, NA),
         col = c(col_a, col_b, col_run), legend = c("ODE", "mean of 2000 runs", "single cells"))
  m10 <- M[4, ]
  n <- 0:25
  plot(n, tabulate(m10 + 1, 26) / 2000, type = "h", lwd = 6, lend = 1, col = adjustcolor(col_a, 0.6),
       xlab = "mRNA molecules at t = 10 h", ylab = "fraction of cells", main = "2000 cells vs Poisson(10)",
       ylim = c(0, 0.15))
  points(n, dpois(n, 10), pch = 19, col = col_b, cex = 1.1)
  legend("topright", bty = "n", cex = 1.2, lwd = c(6, NA), pch = c(NA, 19),
         col = c(adjustcolor(col_a, 0.6), col_b), legend = c("SSA", "Poisson(10)"))
})

## ---- fig05: telegraph trajectories ------------------------------------------

set.seed(5)
pc <- ssa_path(c(1, 10, 100), S3, constit, 40)
pb <- ssa_path(c(0, 10, 100), S3, telegraph, 40)
fig("fig05-telegraph.svg", w = 11, h = 5.6, mfrow = c(2, 2), {
  par(mar = c(4.6, 5.6, 2.4, 1), mgp = c(3, 0.8, 0))
  for (k in 1:2) {
    p <- if (k == 1) pc else pb
    plot(p$t, p$x[, 2], type = "s", col = if (k == 1) col_a else col_b, lwd = 1.4,
         ylim = c(0, 45), xlab = "", ylab = "mRNA",
         main = if (k == 1) "constitutive" else "telegraph (bursty)")
    if (k == 2) {                          # shade the on periods
      on <- p$x[, 1] == 1
      idx <- which(on)
      for (i in idx) rect(p$t[i], -5, p$t[i + 1], 75, col = adjustcolor(col_c, 0.45), border = NA)
      lines(p$t, p$x[, 2], type = "s", col = col_b, lwd = 1.4)
    }
    abline(h = 10, lty = 2, col = col_ref)
  }
  for (k in 1:2) {
    p <- if (k == 1) pc else pb
    plot(p$t, p$x[, 3], type = "s", col = if (k == 1) col_a else col_b, lwd = 1.4,
         ylim = c(0, 250), xlab = "time (h)", ylab = "protein")
    abline(h = 100, lty = 2, col = col_ref)
  }
})

## ---- fig06: protein distributions -------------------------------------------

fig("fig06-protein.svg", w = 11, h = 5.0, mfrow = c(1, 2), {
  br <- seq(0, 500, by = 10)
  for (k in 1:2) {
    v <- if (k == 1) Xc[, 3] else Xb[, 3]
    h <- hist(pmin(v, 499), breaks = br, plot = FALSE)
    plot(h, col = adjustcolor(if (k == 1) col_a else col_b, 0.6), border = "white",
         freq = FALSE, xlab = "protein molecules per cell", ylab = "density", ylim = c(0, 0.042),
         main = if (k == 1) "constitutive" else "telegraph (bursty)")
    curve(dpois(round(x), 100), add = TRUE, n = 501, lwd = 2.5, col = col_ref, lty = 2)
    legend("topright", bty = "n", cex = 1.25,
           legend = c(sprintf("mean %.1f", mean(v)), sprintf("Fano %.1f", fano(v)), "dashed: Poisson(100)"))
  }
})

## ---- fig07: epidemic from one case -------------------------------------------

sir_path <- function(I0, N = 10000, beta = 0.3, gamma = 0.1, t_end = 160) {
  S <- N - I0; I <- I0; t <- 0
  tt <- numeric(25000); II <- numeric(25000); n <- 1; tt[1] <- 0; II[1] <- I
  while (I > 0) {
    a1 <- beta * S * I / N; a2 <- gamma * I; a0 <- a1 + a2
    t <- t - log(runif(1)) / a0
    if (t > t_end) break
    if (runif(1) * a0 < a1) { S <- S - 1; I <- I + 1 } else I <- I - 1
    n <- n + 1; tt[n] <- t; II[n] <- I
  }
  n <- n + 1; tt[n] <- t_end; II[n] <- I
  list(t = tt[1:n], I = II[1:n])
}
sir_ode <- function(I0, h = 0.01) {
  N <- 10000; S <- N - I0; I <- I0; out <- numeric(16001)
  for (n in 1:16001) { out[n] <- I
    dS <- -0.3 * S * I / N; dI <- 0.3 * S * I / N - 0.1 * I
    S <- S + h * dS; I <- I + h * dI }
  list(t = seq(0, 160, by = h), I = out)
}
set.seed(7)
runs <- lapply(1:30, function(r) sir_path(1))
ode <- sir_ode(1)
fig("fig07-sir.png", w = 10, h = 5.2, {
  plot(NA, xlim = c(0, 160), ylim = c(0, 3500), xlab = "time (days)", ylab = "infected I(t)",
       main = "30 simulated outbreaks from I(0) = 1, and the ODE")
  for (r in runs) lines(r$t, r$I, type = "s", col = col_run)
  lines(ode$t, ode$I, lwd = 4, col = col_a)
  died <- sum(sapply(runs, function(r) tail(r$I, 1) == 0 && max(r$I) < 100))
  legend("topright", bty = "n", cex = 1.25, lwd = c(4, 1), col = c(col_a, col_run),
         legend = c("ODE (Euler, h = 0.01 d)", sprintf("SSA runs (%d of 30 die out early)", died)))
})
cat("done\n")
