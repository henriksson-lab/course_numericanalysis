# Figures for the monte-carlo/ unit ("Random numbers and Monte Carlo").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd monte-carlo && Rscript make-monte-carlo-figs.R
#
# Base R graphics only; writes fig01-...svg to fig08-...svg. Everything is
# computed here, with the same code, seeds and numbers as on the slides.

col_a   <- "#0072B2"   # blue
col_b   <- "#D55E00"   # vermillion
col_c   <- "#009E73"   # bluish green
col_ref <- "grey45"
col_pur <- "#CC79A7"   # reddish purple

fig <- function(file, w = 9, h = 5.2, expr, mfrow = c(1, 1)) {
  svg(file, width = w, height = h)
  par(mfrow = mfrow, mar = c(4.8, 5.6, 2.2, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3, cex.main = 1.4)
  expr
  dev.off()
}
logaxis <- function(side, at) {
  axis(side, at = at, labels = parse(text = paste0("10^", log10(at))))
}

## ---- fig01: binding sites, simulated null distribution ----------------------

G <- 1e6; L <- 10
pstart <- seq(10001, by = 20000, length.out = 50)
pend   <- pstart + 999
n_hits <- 200; k_obs <- 18
p <- 50 * (1000 + L - 1) / (G - L + 1)
overlaps <- function(s) {
  i <- findInterval(s + L - 1, pstart)
  i > 0 & pend[pmax(i, 1)] >= s
}
set.seed(1)
K <- replicate(10000, sum(overlaps(sample.int(G - L + 1, n_hits, replace = TRUE))))

fig("fig01-binding-sites.svg", w = 10, h = 5.4, {
  tab <- table(factor(K, levels = 0:30)) / length(K)
  plot(0:30, as.numeric(tab), type = "h", lwd = 7, col = col_a, lend = 1,
       xlab = "number of the 200 hits that overlap a promoter",
       ylab = "probability", ylim = c(0, 0.14))
  points(0:30, dbinom(0:30, n_hits, p), pch = 19, col = col_b, cex = 1.2)
  abline(v = k_obs, lty = 2, lwd = 2.5, col = "black")
  text(k_obs + 0.3, 0.06, "observed: 18", adj = 0, cex = 1.3)
  legend("topright", c("10 000 random placements", "exact binomial"),
         col = c(col_a, col_b), lwd = c(7, NA), pch = c(NA, 19), bty = "n",
         cex = 1.25, inset = c(0, 0.12))
})

## ---- fig02: lattice structure of a bad LCG -------------------------------------

lcg_u <- function(n, seed, a, c, m) {
  x <- numeric(n); x[1] <- seed
  for (i in 2:n) x[i] <- (a * x[i - 1] + c) %% m
  x / m
}
u_bad <- lcg_u(1001, seed = 1, a = 9, c = 1, m = 2048)
set.seed(1); u_good <- runif(1001)
fig("fig02-lcg-lattice.svg", w = 10, h = 5.2, mfrow = c(1, 2), {
  par(mar = c(4.8, 5.2, 3, 1))
  plot(u_bad[-1001], u_bad[-1], pch = 19, cex = 0.5, col = col_b, asp = 1,
       xlab = expression(u[n]), ylab = expression(u[n + 1]),
       main = "LCG, a = 9, m = 2048")
  plot(u_good[-1001], u_good[-1], pch = 19, cex = 0.5, col = col_a, asp = 1,
       xlab = expression(u[n]), ylab = expression(u[n + 1]),
       main = "runif (Mersenne Twister)")
})

## ---- fig03: inverse-transform sampling -----------------------------------------

lambda <- 0.5
fig("fig03-inverse-transform.svg", w = 10, h = 5.2, mfrow = c(1, 2), {
  par(mar = c(4.8, 5.2, 3, 1))
  xs <- seq(0, 10, length.out = 400)
  plot(xs, 1 - exp(-lambda * xs), type = "l", lwd = 3, col = col_a,
       xlab = "waiting time x (h)", ylab = "F(x) = 1 - exp(-0.5 x)",
       main = "draw U, read off x", ylim = c(0, 1))
  us <- c(0.3, 0.75, 0.95)
  for (uu in us) {
    xx <- -log(1 - uu) / lambda
    arrows(0, uu, xx, uu, length = 0.1, col = col_b, lwd = 2)
    arrows(xx, uu, xx, 0, length = 0.1, col = col_b, lwd = 2)
    text(xx + 0.25, uu - 0.05, paste0("U = ", uu), adj = 0, cex = 1.1, col = col_b)
  }
  set.seed(1)
  x <- -log(runif(1e5)) / lambda
  hist(x[x < 12], breaks = seq(0, 12, by = 0.5), freq = FALSE, col = "grey85",
       border = "white", main = "100 000 draws", xlab = "waiting time x (h)",
       ylab = "density")
  curve(lambda * exp(-lambda * x), add = TRUE, lwd = 3, col = col_b)
})

## ---- fig04: Monte Carlo integration as an average height -----------------------

K_ <- 1.2; N0 <- 0.05; r <- 0.8
od <- function(t) K_ / (1 + (K_ / N0 - 1) * exp(-r * t))
exact <- K_ / r * log((exp(r * 10) + K_ / N0 - 1) / (K_ / N0))
fig("fig04-mc-auc.svg", w = 9.5, h = 5.0, {
  curve(od, 0, 10, lwd = 3, col = col_a, xlab = "time (h)", ylab = "OD",
        ylim = c(0, 1.3))
  set.seed(3)
  tt <- runif(20, 0, 10)
  segments(tt, 0, tt, od(tt), col = col_b, lwd = 2)
  points(tt, od(tt), pch = 19, col = col_b)
  abline(h = mean(od(tt)), lty = 2, lwd = 2.5)
  text(0.2, 1.25,
       sprintf("mean height %.3f  x  width 10  =  %.2f  (exact %.2f)",
               mean(od(tt)), 10 * mean(od(tt)), exact), adj = 0, cex = 1.2)
})

## ---- fig05: error against N, Monte Carlo vs trapezoid ---------------------------

mc <- function(N) { t <- runif(N, 0, 10); 10 * mean(od(t)) }
trap <- function(N) { t <- seq(0, 10, length.out = N); h <- t[2] - t[1]
  h * (sum(od(t)) - (od(0) + od(10)) / 2) }
set.seed(1)
Ns <- 10^(2:6)
err_mc <- sapply(Ns, function(N) sqrt(mean((replicate(20, mc(N)) - exact)^2)))
Nt <- 10^(1:4)
err_tr <- abs(sapply(Nt, trap) - exact)
fig("fig05-error-vs-N.svg", w = 9, h = 5.4, {
  plot(Ns, err_mc, log = "xy", type = "b", pch = 19, lwd = 2.5, col = col_a,
       xlim = c(10, 1e6), ylim = c(1e-10, 1), axes = FALSE,
       xlab = "number of function evaluations N", ylab = "error")
  logaxis(1, 10^(1:6)); logaxis(2, 10^seq(-10, 0, by = 2)); box()
  lines(Nt, err_tr, type = "b", pch = 17, lwd = 2.5, col = col_b)
  lines(Ns, 4 * Ns^-0.5, lty = 2, col = col_ref, lwd = 2)
  text(2e3, 0.6, expression("slope " - 1/2), cex = 1.3, col = col_a)
  text(1.5e2, 3e-6, expression("trapezoid: slope " - 2), cex = 1.3, col = col_b, adj = 0)
  legend("bottomleft", c("Monte Carlo (RMS of 20 runs)", "trapezoid rule"),
         col = c(col_a, col_b), pch = c(19, 17), lwd = 2.5, bty = "n", cex = 1.2)
})

## ---- fig06: running estimates, finite vs infinite variance ----------------------

run_est <- function(fx) {
  n <- seq_along(fx)
  m <- cumsum(fx) / n
  s <- sqrt(pmax(cumsum(fx^2) / n - m^2, 0) * n / pmax(n - 1, 1))
  list(n = n, m = m, se = s / sqrt(n))
}
set.seed(1)
x <- runif(1e5)
ra <- run_est(10 * od(10 * x))
rb <- run_est(x^(-0.75))
fig("fig06-heavy-tail.svg", w = 10, h = 5.2, mfrow = c(1, 2), {
  par(mar = c(4.8, 5.2, 3, 1))
  ii <- unique(round(10^seq(1, 5, length.out = 600)))
  plot(ra$n[ii], ra$m[ii], log = "x", type = "l", lwd = 2.5, col = col_a,
       ylim = c(6.5, 8), xlab = "N", ylab = "running estimate", axes = FALSE,
       main = "OD curve: finite variance")
  logaxis(1, 10^(1:5)); axis(2); box()
  polygon(c(ra$n[ii], rev(ra$n[ii])),
          c(ra$m[ii] - 2 * ra$se[ii], rev(ra$m[ii] + 2 * ra$se[ii])),
          col = adjustcolor(col_a, 0.2), border = NA)
  abline(h = exact, lty = 2, lwd = 2)
  plot(rb$n[ii], rb$m[ii], log = "x", type = "l", lwd = 2.5, col = col_b,
       ylim = c(2, 6), xlab = "N", ylab = "running estimate", axes = FALSE,
       main = expression(x^{-3/4} * ": infinite variance"))
  logaxis(1, 10^(1:5)); axis(2); box()
  polygon(c(rb$n[ii], rev(rb$n[ii])),
          c(rb$m[ii] - 2 * rb$se[ii], rev(rb$m[ii] + 2 * rb$se[ii])),
          col = adjustcolor(col_b, 0.2), border = NA)
  abline(h = 4, lty = 2, lwd = 2)
})

## ---- fig07: curse of dimensionality --------------------------------------------
# Toy marginal likelihood: d parameters, prior N(0,1) each, one measurement
# y_i = 1 with noise sd 1. Z = dnorm(1, 0, sqrt(2))^d exactly. The grid is the
# midpoint rule in u = pnorm(theta) with n points per axis, n^d <= 1e5. Because
# the integrand is a product, the grid sum equals the product of 1-D sums,
# which is what we compute here (the slides build the full grid for d = 10).

Z1 <- dnorm(1, 0, sqrt(2))
ds <- 1:20
res <- t(sapply(ds, function(d) {
  n <- floor(1e5^(1 / d) + 1e-9); N <- n^d
  nodes <- qnorm((1:n - 0.5) / n)
  Zg <- mean(dnorm(1, nodes, 1))^d
  set.seed(1)
  th <- matrix(rnorm(d * 1e5), nrow = d)
  Zm <- mean(exp(colSums(dnorm(1, th, 1, log = TRUE))))
  c(d = d, n = n, grid = abs(Zg / Z1^d - 1), mc = abs(Zm / Z1^d - 1))
}))
fig("fig07-curse.svg", w = 9.5, h = 5.4, {
  plot(res[, "d"], res[, "grid"], log = "y", type = "b", pch = 17, lwd = 2.5,
       col = col_b, ylim = c(1e-10, 1e3), axes = FALSE,
       xlab = "dimension d", ylab = "relative error")
  logaxis(2, 10^seq(-10, 2, by = 2)); axis(1, at = c(1, 5, 10, 15, 20)); box()
  lines(res[, "d"], res[, "mc"], type = "b", pch = 19, lwd = 2.5, col = col_a)
  sel <- c(1, 2, 3, 5, 10, 17)
  text(res[sel, "d"], res[sel, "grid"], paste0("n=", format(res[sel, "n"], scientific = FALSE, trim = TRUE)),
       pos = 4, cex = 1.05, col = col_b)
  legend("bottomright", c("grid, n points per axis, n^d <= 100 000",
                          "Monte Carlo, 100 000 draws"),
         col = c(col_b, col_a), pch = c(17, 19), lwd = 2.5, bty = "n", cex = 1.2)
})

## ---- fig08: rare event and the importance-sampling proposal --------------------

fig("fig08-importance.svg", w = 9.5, h = 5.0, {
  xs <- seq(-4, 9, length.out = 600)
  plot(xs, dnorm(xs), type = "l", lwd = 3, col = col_a, log = "y",
       ylim = c(1e-12, 1), axes = FALSE, xlab = "standardised motif score x",
       ylab = "density")
  logaxis(2, 10^seq(-12, 0, by = 3)); axis(1); box()
  lines(xs, dnorm(xs, mean = 5), lwd = 3, col = col_b, lty = 2)
  xt <- xs[xs >= 5 & dnorm(xs) > 1e-12]
  polygon(c(5, xt, max(xt)), c(1e-12, dnorm(xt), 1e-12),
          col = adjustcolor(col_a, 0.3), border = NA)
  abline(v = 5, lwd = 1.5, col = col_ref)
  mtext("f = N(0, 1): where plain MC draws", side = 3, line = 0.5, adj = 0,
        cex = 1.2, col = col_a)
  mtext("g = N(5, 1): proposal (dashed)", side = 3, line = 0.5, adj = 1,
        cex = 1.2, col = col_b)
  text(5.6, 3e-9, "P(X > 5) = 2.9e-7", adj = 0, cex = 1.2)
})
