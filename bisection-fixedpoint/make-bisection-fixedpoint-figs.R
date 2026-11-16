# Figures for the bisection-fixedpoint/ unit ("Bisection and fixed-point iteration").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd bisection-fixedpoint && Rscript make-bisection-fixedpoint-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Makes:
#   fig01-harvest-roots.svg   f(N) = g(N) - N for the harvested Ricker stock, two roots
#   fig02-no-root.svg         f(N) for catches H = 10, 14.26 (maximum) and 16: roots merge and vanish
#   fig03-cobweb-bh.svg       cobweb plot of the Beverton-Holt model (R = 2, K = 100)
#   fig04-cobweb-harvest.svg  cobweb of the harvested Ricker stock from N0 = 21 and N0 = 20
#   fig05-error-iter.svg      error vs iteration: bisection, population map g, rearranged g2
#   fig06-ricker-regimes.svg  cobwebs and time series of the Ricker model for r = 0.5, 1.5, 2.3, 3
#   fig07-bifurcation.svg     Ricker bifurcation diagram: long-run values against r

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- the running example (same code as on the slides) ----------------------

r <- 0.5; K <- 100; H <- 10            # growth rate, capacity, catch (thousand fish)
g  <- function(N, H = 10) N * exp(r * (1 - N / K)) - H   # stock next year
f  <- function(N, H = 10) g(N, H) - N                    # steady state: f(N) = 0
g2 <- function(N) K * (1 - log(1 + H / N) / r)           # rearranged: converges to upper root
root_lo <- uniroot(f, c(0, 50),   tol = 1e-14)$root      # 20.48199
root_hi <- uniroot(f, c(50, 100), tol = 1e-14)$root      # 74.95244

setup <- function(mar = c(4.8, 5.4, 1.5, 1)) {
  par(mar = mar, las = 1, mgp = c(3.4, 0.8, 0), cex.lab = 1.5, cex.axis = 1.35)
}

# cobweb: vertical to the curve, horizontal to the diagonal
cobweb <- function(gfun, x0, n, col, lwd = 2) {
  x <- x0
  for (k in 1:n) {
    y <- gfun(x)
    segments(x, if (k == 1) 0 else x, x, y, col = col, lwd = lwd)
    segments(x, y, y, y, col = col, lwd = lwd)
    x <- y
  }
}

## ---- fig01: two roots of f -------------------------------------------------

svg("fig01-harvest-roots.svg", width = 9, height = 5.2)
setup()
NN <- seq(0, 100, length.out = 501)
plot(NN, f(NN), type = "n", xlab = "stock N (thousand fish)",
     ylab = "f(N) = g(N) - N", ylim = c(-10.5, 5.5))
rect(50, -11.5, 100, 6.5, col = adjustcolor(oi["sky"], 0.15), border = NA)
abline(h = 0, col = oi["grey"], lwd = 1.5)
lines(NN, f(NN), col = oi["blue"], lwd = 3)
points(c(root_lo, root_hi), c(0, 0), pch = 19, cex = 1.8, col = oi["red"])
text(root_lo, 1.2, sprintf("%.2f", root_lo), cex = 1.35, col = oi["red"])
text(root_hi, 1.2, sprintf("%.2f", root_hi), cex = 1.35, col = oi["red"], pos = 4)
points(c(50, 100), f(c(50, 100)), pch = 15, cex = 1.6, col = oi["orange"])
text(50, f(50), "f(50) > 0", pos = 3, cex = 1.3, col = oi["orange"])
text(100, f(100), "f(100) < 0", pos = 2, cex = 1.3, col = oi["orange"])
text(75, -8, "bracket [50, 100]", cex = 1.35, col = oi["blue"])
text(2, 4.6, "f > 0: stock grows", cex = 1.2, col = oi["grey"], adj = 0)
text(2, -2.5, "f < 0:\nstock shrinks", cex = 1.2, col = oi["grey"], adj = 0)
dev.off()

## ---- fig02: increasing the catch: roots merge and vanish --------------------

svg("fig02-no-root.svg", width = 9, height = 5.2)
setup()
s <- function(N) N * (exp(r * (1 - N / K)) - 1)          # surplus production
Hmax <- optimize(s, c(0, 100), maximum = TRUE)$objective
Hs <- c(10, Hmax, 16)
cols <- c(oi["blue"], oi["orange"], oi["red"])
plot(NN, f(NN, 10), type = "n", xlab = "stock N (thousand fish)",
     ylab = "f(N) = g(N) - N", ylim = c(-22, 5))
abline(h = 0, col = oi["grey"], lwd = 1.5)
for (i in 1:3) lines(NN, f(NN, Hs[i]), col = cols[i], lwd = 3, lty = c(1, 2, 1)[i])
legend("bottom", lwd = 3, lty = c(1, 2, 1), col = cols, bty = "n", cex = 1.3, horiz = TRUE,
       legend = c("H = 10: two roots", sprintf("H = %.2f: one (double)", Hmax),
                  "H = 16: none"))
dev.off()

## ---- fig03: Beverton-Holt cobweb --------------------------------------------

svg("fig03-cobweb-bh.svg", width = 7, height = 6)
setup()
R0 <- 2
bh <- function(N) R0 * N / (1 + (R0 - 1) * N / K)
x <- seq(0, 130, length.out = 300)
plot(x, bh(x), type = "n", xlab = expression(N[t]), ylab = expression(N[t + 1]),
     xlim = c(0, 130), ylim = c(0, 130), xaxs = "i", yaxs = "i")
abline(0, 1, col = oi["grey"], lwd = 2)
lines(x, bh(x), col = oi["blue"], lwd = 3)
cobweb(bh, 20, 8, oi["red"])
points(K, K, pch = 19, cex = 1.8)
text(K, K, "  N* = K", pos = 4, cex = 1.3)
text(20, 4, expression(N[0]), cex = 1.3, pos = 4)
text(75, 62, expression("diagonal" ~ N[t + 1] == N[t]), cex = 1.2, col = oi["grey"], adj = 0)
dev.off()

## ---- fig04: cobweb of the harvested stock, two starts -----------------------

svg("fig04-cobweb-harvest.svg", width = 7.5, height = 6)
setup()
x <- seq(0, 110, length.out = 400)
plot(x, g(x), type = "n", xlab = expression(N[t]), ylab = expression(N[t + 1]),
     xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i")
abline(0, 1, col = oi["grey"], lwd = 2)
lines(x, g(x), col = oi["blue"], lwd = 3)
cobweb(g, 21, 40, oi["green"])
gpos <- function(N) pmax(g(N), 0)                        # stop drawing at extinction
cobweb(gpos, 20, 12, oi["red"])
points(c(root_lo, root_hi), c(root_lo, root_hi), pch = c(1, 19), cex = 2, lwd = 2)
text(root_hi, root_hi - 6, "stable\ng' = 0.71", cex = 1.2, pos = 4)
text(30, 14, "unstable\ng' = 1.34", cex = 1.2, adj = 0)
arrows(29, 16, root_lo + 1.5, root_lo - 1, length = 0.1)
legend("topleft", lwd = 2.5, col = c(oi["green"], oi["red"]), bty = "n", cex = 1.25,
       legend = c(expression(N[0] == 21 ~ "→ 74.95"), expression(N[0] == 20 ~ "→ collapse")))
dev.off()

## ---- fig05: error vs iteration ----------------------------------------------

svg("fig05-error-iter.svg", width = 9, height = 5.4)
setup()
n <- 60
# bisection on [50, 100]: error of the midpoint and the bound (b - a)/2
a <- 50; b <- 100; eb <- bnd <- numeric(n)
for (k in 1:n) {
  m <- (a + b) / 2
  eb[k] <- abs(m - root_hi); bnd[k] <- (b - a) / 2
  if (sign(f(m)) == sign(f(a))) a <- m else b <- m
}
iter_err <- function(gfun, x0) {
  e <- numeric(n); x <- x0
  for (k in 1:n) { x <- gfun(x); e[k] <- abs(x - root_hi) }
  e
}
eg  <- iter_err(g, 30)
eg2 <- iter_err(g2, 30)
clip <- function(e) ifelse(e < 1e-16, NA, e)
plot(1:n, clip(bnd), log = "y", type = "n", ylim = c(1e-16, 1e5), yaxt = "n",
     xlab = "iteration k", ylab = "error")
axis(2, at = 10^seq(-15, 0, by = 3),
     labels = parse(text = paste0("10^", seq(-15, 0, by = 3))))
lines(1:n, clip(bnd), col = oi["orange"], lwd = 3)
points(1:n, clip(eb), col = oi["orange"], pch = 1, cex = 0.8)
lines(1:n, clip(eg), col = oi["blue"], lwd = 3)
lines(1:n, clip(eg2), col = oi["green"], lwd = 3)
legend("topright", bty = "n", cex = 1.15, lwd = c(3, NA, 3, 3), pch = c(NA, 1, NA, NA),
       col = c(oi["orange"], oi["orange"], oi["blue"], oi["green"]),
       legend = c("bisection: bound (b - a)/2, factor 0.5",
                  "bisection: actual midpoint error",
                  "population map g: factor |g'| = 0.71",
                  "rearranged g2: factor |g2'| = 0.31"))
dev.off()

## ---- fig06: Ricker regimes ---------------------------------------------------

ricker <- function(N, rr) N * exp(rr * (1 - N / K))
svg("fig06-ricker-regimes.svg", width = 13, height = 6.2)
par(mfrow = c(2, 4), mar = c(4.2, 4.6, 3, 0.8), las = 1, mgp = c(2.8, 0.7, 0),
    cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.6)
rs <- c(0.5, 1.5, 2.3, 3.0)
lab <- c("monotone: g' = 0.5", "oscillating: g' = -0.5",
         "2-cycle: g' = -1.3", "chaos: g' = -2")
x <- seq(0, 260, length.out = 500)
for (i in 1:4) {
  rr <- rs[i]
  plot(x, ricker(x, rr), type = "n", xlab = expression(N[t]), ylab = expression(N[t + 1]),
       xlim = c(0, 250), ylim = c(0, 250), xaxs = "i", yaxs = "i",
       main = paste0("r = ", rr))
  abline(0, 1, col = oi["grey"], lwd = 1.5)
  lines(x, ricker(x, rr), col = oi["blue"], lwd = 2.5)
  cobweb(function(N) ricker(N, rr), 20, 40, oi["red"], lwd = 1.3)
  points(K, K, pch = 19, cex = 1.4)
}
for (i in 1:4) {
  rr <- rs[i]
  N <- numeric(41); N[1] <- 20
  for (t in 1:40) N[t + 1] <- ricker(N[t], rr)
  plot(0:40, N, type = "o", pch = 19, cex = 0.6, col = oi["red"], ylim = c(0, 250),
       xlab = "year t", ylab = expression(N[t]))
  abline(h = K, col = oi["grey"], lty = 2)
  mtext(lab[i], side = 3, line = 0.4, cex = 1.0)
}
dev.off()

## ---- fig07: bifurcation diagram ----------------------------------------------

svg("fig07-bifurcation.svg", width = 10, height = 5.4)
setup()
rgrid <- seq(0.05, 3.4, length.out = 1200)
plot(NA, xlim = range(rgrid), ylim = c(0, 420), xlab = "growth rate r",
     ylab = "long-run stock N")
rect(0, -50, 2, 500, col = adjustcolor(oi["green"], 0.12), border = NA)
for (rr in rgrid) {
  N <- 20
  for (t in 1:500) N <- ricker(N, rr)
  keep <- numeric(80)
  for (t in 1:80) { N <- ricker(N, rr); keep[t] <- N }
  points(rep(rr, 80), keep, pch = ".", col = oi["blue"], cex = 1.6)
}
segments(2, K, 3.4, K, col = oi["red"], lty = 2, lwd = 2)
text(1, 380, "|g'(K)| = |1 - r| < 1\nsteady state stable", cex = 1.3)
text(2.75, 380, "|1 - r| > 1:\ncycles, then chaos", cex = 1.3)
text(2.25, K + 14, "K (unstable)", col = oi["red"], cex = 1.15, adj = 0)
dev.off()
