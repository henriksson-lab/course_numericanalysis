# Figures for the nonlinear-least-squares/ unit
# ("Nonlinear least squares and ODE fitting").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd nonlinear-least-squares && Rscript make-nonlinear-least-squares-figs.R
#
# Base R graphics plus deSolve and minpack.lm; writes SVG files fig01-...svg to
# fig08-...svg. The data are the simulated plasma concentrations used on the
# slides and in the exercise: one-compartment model with first-order oral
# absorption, dose 100 mg, ka = 1 /h, ke = 0.2 /h, V = 20 L, additive noise
# N(0, 0.1^2) mg/L (set.seed(1)), rounded to 2 decimals and typed in.

suppressMessages({library(deSolve); library(minpack.lm)})

## ---- the running example (same numbers as on the slides) ---------------------

D  <- 100
tt <- c(0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8, 10, 12)
y  <- c(1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61)
bateman <- function(t, ka, ke, V) D * ka / (V * (ka - ke)) * (exp(-ke * t) - exp(-ka * t))
model <- function(p) bateman(tt, p[1], p[2], p[3])
rss <- function(p) sum((y - model(p))^2)

jac <- function(p, h = 1e-7) sapply(1:3, function(j) {
  e <- replace(numeric(3), j, h * p[j]); (model(p + e) - model(p)) / e[j] })

gn_path <- function(p, maxit = 30, tol = 1e-8) {
  path <- rbind(p)
  for (k in 1:maxit) {
    delta <- tryCatch(qr.solve(jac(p), y - model(p)), error = function(e) NULL)
    if (is.null(delta)) break
    p <- p + delta; path <- rbind(path, p)
    if (!all(is.finite(model(p))) || sqrt(sum(delta^2)) < tol * sqrt(sum(p^2))) break
  }
  path
}
lm_path <- function(p, lambda = 0.01, maxit = 200, tol = 1e-8) {
  path <- rbind(p); rss0 <- rss(p)
  for (k in 1:maxit) {
    J <- jac(p); r <- y - model(p)
    delta <- qr.solve(rbind(J, sqrt(lambda) * diag(sqrt(colSums(J^2)))), c(r, 0, 0, 0))
    if (sqrt(sum(delta^2)) < tol * sqrt(sum(p^2))) break
    new <- rss(p + delta)
    if (is.finite(new) && new < rss0) {
      p <- p + delta; rss0 <- new; lambda <- lambda / 10; path <- rbind(path, p)
    } else lambda <- lambda * 10
  }
  path
}

p_hat <- coef(nlsLM(y ~ bateman(tt, ka, ke, V), start = list(ka = 0.5, ke = 0.1, V = 10),
                    control = nls.lm.control(ftol = 1e-15, ptol = 1e-15)))
p_flip <- c(p_hat[2], p_hat[1], p_hat[3] * p_hat[2] / p_hat[1])

# RSS with the scale D/V chosen optimally for each (ka, ke) ("profiled")
rss_prof <- function(ka, ke) {
  g <- ka / (ka - ke) * (exp(-ke * tt) - exp(-ka * tt))
  a <- sum(g * y) / sum(g^2)
  sum((y - a * g)^2)
}

## ---- plotting helpers ---------------------------------------------------------

col_a <- "#0072B2"   # blue
col_b <- "#D55E00"   # vermillion
col_c <- "#009E73"   # bluish green
col_d <- "#E69F00"   # orange
col_ref <- "grey45"

fig <- function(file, w = 9, h = 5.2, expr) {
  svg(file, width = w, height = h)
  par(mar = c(4.8, 5.6, 1.5, 1), las = 1, mgp = c(3.6, 0.8, 0),
      cex.lab = 1.4, cex.axis = 1.3)
  expr
  dev.off()
}
logaxis <- function(side, at) {
  axis(side, at = 10^at, labels = parse(text = paste0("10^", at)))
}

## ---- fig01: the data and the fitted curve ---------------------------------------

fig("fig01-pk-data.svg", expr = {
  tg <- seq(0, 12.5, length.out = 400)
  plot(tt, y, pch = 19, col = col_a, cex = 1.4, xlim = c(0, 12.5), ylim = c(0, 3.8),
       xlab = "time after oral dose (h)", ylab = "plasma concentration (mg/L)")
  lines(tg, bateman(tg, 1, 0.2, 20), col = col_ref, lwd = 2, lty = 2)
  lines(tg, bateman(tg, p_hat[1], p_hat[2], p_hat[3]), col = col_b, lwd = 2.5)
  legend("topright", c("measured (simulated)", "true: ka = 1, ke = 0.2, V = 20",
                       "least-squares fit"),
         pch = c(19, NA, NA), lty = c(NA, 2, 1), lwd = c(NA, 2, 2.5),
         col = c(col_a, col_ref, col_b), bty = "n", cex = 1.25)
})

## ---- fig02: the RSS landscape over (ka, ke), V profiled out -------------------

lg <- seq(log10(0.03), log10(10), length.out = 241)
Z <- outer(10^lg, 10^lg, Vectorize(function(a, e) if (abs(a - e) < 1e-9) NA else rss_prof(a, e)))
fig("fig02-rss-landscape.svg", w = 7.4, h = 6, expr = {
  par(mar = c(4.8, 5.8, 1.5, 1))
  lev <- c(0.07, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40)
  image(lg, lg, log10(Z), col = hcl.colors(60, "Blues 3", rev = TRUE), axes = FALSE, useRaster = TRUE,
        xlab = expression(k[a] ~ "(per hour)"), ylab = expression(k[e] ~ "(per hour)"))
  contour(lg, lg, Z, levels = lev, add = TRUE, labcex = 0.9, col = "grey25")
  ax <- c(0.03, 0.1, 0.3, 1, 3, 10)
  axis(1, at = log10(ax), labels = ax); axis(2, at = log10(ax), labels = ax); box()
  abline(0, 1, lty = 2, lwd = 2, col = col_ref)
  points(log10(p_hat[1]), log10(p_hat[2]), pch = 8, cex = 2.2, lwd = 2.5, col = col_b)
  points(log10(p_flip[1]), log10(p_flip[2]), pch = 8, cex = 2.2, lwd = 2.5, col = col_d)
  text(log10(p_hat[1]), log10(p_hat[2]) - 0.17, "V = 20.2 L", cex = 1.2, col = col_b, font = 2)
  text(log10(p_flip[1]), log10(p_flip[2]) + 0.17, "V = 3.8 L", cex = 1.2, col = col_d, font = 2)
  text(log10(5), log10(2.2), expression(k[a] == k[e]), cex = 1.3, col = "grey20")
})

## ---- fig03: convergence of Gauss-Newton and Levenberg-Marquardt ---------------

relerr <- function(path) apply(path, 1, function(p) sqrt(sum((p - p_hat)^2)) / sqrt(sum(p_hat^2)))
e_gn  <- relerr(gn_path(c(0.5, 0.1, 10)))
e_lm  <- relerr(lm_path(c(0.5, 0.1, 10)))
e_lm2 <- relerr(lm_path(c(3, 0.05, 50)))
fig("fig03-convergence.svg", expr = {
  plot(seq_along(e_gn) - 1, e_gn, log = "y", type = "o", pch = 19, lwd = 2.5, col = col_a,
       xlim = c(0, 12), ylim = c(1e-12, 3), yaxt = "n",
       xlab = "accepted iteration", ylab = "relative error in (ka, ke, V)")
  logaxis(2, seq(-12, 0, by = 3))
  lines(seq_along(e_lm) - 1, e_lm, type = "o", pch = 17, lwd = 2.5, col = col_c)
  lines(seq_along(e_lm2) - 1, e_lm2, type = "o", pch = 15, lwd = 2.5, col = col_b)
  legend("bottomleft", c("Gauss-Newton, start (0.5, 0.1, 10)",
                         "Levenberg-Marquardt, same start",
                         "Levenberg-Marquardt, start (3, 0.05, 50)"),
         col = c(col_a, col_c, col_b), pch = c(19, 17, 15), lwd = 2.5, bty = "n", cex = 1.2)
})

## ---- fig04: GN and LM paths from a bad start ------------------------------------

xa <- seq(0.05, 3.5, length.out = 200); ye <- seq(0.01, 0.45, length.out = 200)
Z2 <- outer(xa, ye, Vectorize(function(a, e) if (abs(a - e) < 1e-6) NA else rss_prof(a, e)))
P_lm <- lm_path(c(3, 0.05, 50)); P_gn <- gn_path(c(3, 0.05, 50))
fig("fig04-gn-lm-paths.svg", expr = {
  contour(xa, ye, Z2, levels = c(0.07, 0.1, 0.2, 0.5, 1, 2, 4, 8), labcex = 0.9,
          col = "grey55", xlab = expression(k[a] ~ "(per hour)"),
          ylab = expression(k[e] ~ "(per hour)"), xlim = c(0, 3.5), ylim = c(0, 0.45))
  lines(P_lm[, 1], P_lm[, 2], type = "o", pch = 19, lwd = 2.5, col = col_c)
  arrows(P_gn[1, 1], P_gn[1, 2], 0.02, 0.05 + (P_gn[2, 2] - 0.05) * (3 - 0.02) / (3 - P_gn[2, 1]),
         lwd = 2.5, col = col_b, length = 0.15)
  text(0.05, 0.025, "GN step to ka = -2.2, V = 0.04", col = col_b, cex = 1.15, adj = 0, font = 2)
  points(3, 0.05, pch = 15, cex = 1.8)
  text(3, 0.05, "start", pos = 1, cex = 1.2)
  points(p_hat[1], p_hat[2], pch = 8, cex = 2, lwd = 2.5, col = col_a)
  legend("topright", c("Levenberg-Marquardt", "Gauss-Newton (first step)"),
         col = c(col_c, col_b), lwd = 2.5, pch = c(19, NA), bty = "n", cex = 1.2)
})

## ---- fig05: the RSS seen through an ODE solver -----------------------------------

pk_rhs <- function(t, A, p) list(c(-p[1] * A[1], p[1] * A[1] - p[2] * A[2]))
pk_ode <- function(p, tol) {
  out <- ode(c(D, 0), c(0, tt), pk_rhs, p, rtol = tol, atol = tol)
  out[-1, 3] / p[3]
}
kag <- seq(0.9, 1.1, length.out = 801)
dev_tol <- function(tol) sapply(kag, function(k) {
  p <- c(k, p_hat[2:3]); sum((y - pk_ode(p, tol))^2) - rss(p) })
d3 <- dev_tol(1e-3); d6 <- dev_tol(1e-6)
fig("fig05-solver-error.svg", w = 11, h = 5, expr = {
  par(mfrow = c(1, 2), mar = c(4.8, 6.6, 2.6, 1), mgp = c(4.6, 0.8, 0))
  plot(kag, d3 * 1e5, type = "l", lwd = 2.5, col = col_b, xlab = expression(k[a]),
       ylab = expression("RSS"[ode] - "RSS"[exact] ~ (10^-5)), main = "rtol = atol = 1e-3",
       cex.main = 1.4)
  abline(v = p_hat[1], lty = 3, col = col_ref)
  plot(kag, d6 * 1e7, type = "l", lwd = 2.5, col = col_a, xlab = expression(k[a]),
       ylab = expression("RSS"[ode] - "RSS"[exact] ~ (10^-7)), main = "rtol = atol = 1e-6",
       cex.main = 1.4)
  abline(v = p_hat[1], lty = 3, col = col_ref)
})

## ---- fig06: inner tolerance vs finite-difference step ----------------------------

resid_ode <- function(q, tol) y - pk_ode(exp(q), tol)
tols <- 10^-(3:12)
epsf <- c(0, 1e-10)        # FD step sqrt(epsfcn) * |q|; 0 means sqrt(machine eps)
E <- sapply(epsf, function(ef) sapply(tols, function(tol) {
  f <- nls.lm(log(c(0.5, 0.1, 10)), fn = resid_ode, tol = tol,
              control = nls.lm.control(epsfcn = ef))
  abs(exp(f$par[1]) / p_hat[1] - 1)
}))
fig("fig06-tolerance-step.svg", expr = {
  matplot(tols, E, log = "xy", type = "o", pch = c(19, 17), lty = 1, lwd = 2.5,
          col = c(col_b, col_a), xaxt = "n", yaxt = "n", ylim = c(1e-8, 1),
          xlim = c(1e-3, 1e-12),
          xlab = "ODE solver tolerance (rtol = atol)", ylab = expression("relative error in" ~ k[a]))
  logaxis(1, -(3:12)); logaxis(2, seq(-8, 0, by = 2))
  legend("topright", c("FD step 1.5e-8 |q|  (MINPACK default)", "FD step 1e-5 |q|  (epsfcn = 1e-10)"),
         col = c(col_b, col_a), pch = c(19, 17), lwd = 2.5, bty = "n", cex = 1.2)
})

## ---- fig07: multistart -------------------------------------------------------------

resid_log <- function(q) y - bateman(tt, exp(q[1]), exp(q[2]), exp(q[3]))
g <- c(0.03, 0.1, 0.3, 1, 3, 10)
starts <- expand.grid(ka = g, ke = 1.5 * g)
ends <- t(apply(starts, 1, function(s) {
  f <- nls.lm(log(c(s, 10)), fn = resid_log, control = nls.lm.control(maxiter = 200))
  c(exp(f$par), f$deviance)
}))
cls <- ifelse(ends[, 4] > 0.1, 3, ifelse(abs(ends[, 1] - p_hat[1]) < 1e-3, 1, 2))
fig("fig07-multistart.svg", w = 7.4, h = 6, expr = {
  par(mar = c(4.8, 5.8, 1.5, 1))
  plot(log10(starts$ka), log10(starts$ke), pch = c(19, 17, 4)[cls], cex = 1.8, lwd = 2.5,
       col = c(col_b, col_d, "grey30")[cls], axes = FALSE, xlim = log10(c(0.02, 15)),
       ylim = log10(c(0.03, 120)), xlab = expression("start" ~ k[a]),
       ylab = expression("start" ~ k[e]))
  ax <- c(0.03, 0.1, 0.3, 1, 3, 10)
  axis(1, at = log10(ax), labels = ax); axis(2, at = log10(ax), labels = ax); box()
  abline(0, 0, lty = 0); abline(a = 0, b = 1, lty = 2, lwd = 2, col = col_ref)
  legend("topleft", c(sprintf("-> ka = 1.02, ke = 0.19, V = 20.2  (%d)", sum(cls == 1)),
                      sprintf("-> ka = 0.19, ke = 1.02, V = 3.8  (%d)", sum(cls == 2)),
                      sprintf("-> runs off, RSS = 64.5  (%d)", sum(cls == 3))),
         pch = c(19, 17, 4), col = c(col_b, col_d, "grey30"), pt.lwd = 2.5, bty = "o",
         bg = "white", box.col = "white", cex = 1.1)
})

## ---- fig08: the flip-flop: same plasma curve, different gut ----------------------

tg <- seq(0, 12, length.out = 400)
fig("fig08-flipflop.svg", w = 11, h = 5, expr = {
  par(mfrow = c(1, 2), mar = c(4.8, 5.6, 2.6, 1))
  plot(tt, y, pch = 19, col = "grey40", cex = 1.3, ylim = c(0, 3.8),
       xlab = "time (h)", ylab = "concentration (mg/L)", main = "plasma (measured)",
       cex.main = 1.4)
  lines(tg, bateman(tg, p_hat[1], p_hat[2], p_hat[3]), col = col_b, lwd = 4)
  lines(tg, bateman(tg, p_flip[1], p_flip[2], p_flip[3]), col = col_d, lwd = 2.5, lty = 2)
  legend("topright", c("ka = 1.02, ke = 0.19, V = 20.2", "ka = 0.19, ke = 1.02, V = 3.8"),
         col = c(col_b, col_d), lwd = c(4, 2.5), lty = c(1, 2), bty = "n", cex = 1.1)
  plot(tg, D * exp(-p_hat[1] * tg), type = "l", col = col_b, lwd = 3, ylim = c(0, 100),
       xlab = "time (h)", ylab = "amount in gut (mg)", main = "gut (not measured)",
       cex.main = 1.4)
  lines(tg, D * exp(-p_flip[1] * tg), col = col_d, lwd = 3, lty = 2)
})

cat("p_hat:", format(p_hat, digits = 10), "\n")
cat("flip:", format(p_flip, digits = 10), "\n")
cat("rss:", rss(p_hat), rss(p_flip), "\n")
print(signif(E, 2))
print(table(cls))
