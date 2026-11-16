# Figures for the optimisation-1d/ unit ("Optimisation in one dimension").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd optimisation-1d && Rscript make-optimisation-1d-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Makes:
#   fig01-min-vs-root.svg   minus log-likelihood (minimum) and score (root), ZTP counts
#   fig02-decay-data.svg    mRNA decay time course and the fitted exponential
#   fig03-decay-sse.svg     SSE(k) for the decay rate: minimum, convex and concave parts
#   fig04-golden.svg        three golden-section steps on SSE(k)
#   fig05-convergence.svg   error vs function evaluations: golden, Brent, Newton
#   fig06-parabola.svg      parabolic interpolation: a good and a bad three-point start
#   fig07-flat-bottom.svg   the flat bottom: -loglik vs score near the MLE, at 1e-8 scale
#   fig08-newton-min.svg    Newton's local parabola: convex (good) and concave (bad) start
#   fig09-circadian.svg     circadian period: data and RSS(T) with two local minima

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- the running examples (same code as on the slides) ---------------------

x   <- rep(1:5, times = c(52, 30, 12, 5, 1))   # UMI counts, cells with >= 1 UMI
n   <- length(x); sx <- sum(x)
S   <- function(l) sx / l - n / (-expm1(-l))   # score
nll <- function(l) -(sx * log(l) - n * l - n * log(-expm1(-l)))  # minus log-likelihood
lhat <- 1.2184695258301268

t <- c(0, 0.5, 1, 2, 4, 6, 8)                            # hours after transcription stop
y <- c(1.000, 0.804, 0.684, 0.484, 0.217, 0.094, 0.083)  # fraction of mRNA remaining
sse <- function(k) sum((y - exp(-k * t))^2)
d1  <- function(k) sum(2 * (y - exp(-k * t)) * t * exp(-k * t))        # SSE'
d2  <- function(k) sum(2 * t^2 * exp(-k * t) * (2 * exp(-k * t) - y))  # SSE''
khat <- uniroot(d1, c(0.2, 0.6), tol = 1e-15)$root                       # 0.3753873094
kinfl <- uniroot(d2, c(0.4, 2), tol = 1e-12)$root                      # 0.802
ssev <- function(k) sapply(k, sse)

gr <- (sqrt(5) - 1) / 2

## ---- fig01: minimum of -loglik = root of the score ----------------------------

svg("fig01-min-vs-root.svg", width = 12, height = 4.8)
par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
ll <- seq(0.4, 3, length.out = 300)
plot(ll, sapply(ll, nll), type = "l", lwd = 3, col = oi["blue"],
     xlab = expression("rate " * lambda), ylab = expression(-"log-likelihood"),
     main = "minimise  f = -log L")
abline(v = lhat, lty = 3, lwd = 2)
points(lhat, nll(lhat), pch = 19, cex = 1.6, col = oi["red"])
text(lhat + 0.08, nll(lhat) + 4, expression(hat(lambda) == 1.218), adj = 0, cex = 1.4, col = oi["red"])
plot(ll, S(ll), type = "l", lwd = 3, col = oi["orange"],
     xlab = expression("rate " * lambda), ylab = expression("score " * S(lambda) == -f * "'"),
     main = "or find the root of  f '")
abline(h = 0, col = oi["grey"], lwd = 1.5); abline(v = lhat, lty = 3, lwd = 2)
points(lhat, 0, pch = 19, cex = 1.6, col = oi["red"])
dev.off()

## ---- fig02: decay data and fit ---------------------------------------------

svg("fig02-decay-data.svg", width = 9, height = 5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
tt <- seq(0, 8.5, length.out = 300)
plot(t, y, pch = 19, cex = 1.6, col = oi["blue"], ylim = c(0, 1.05), xlim = c(0, 8.5),
     xlab = "time after transcription stop (h)", ylab = "fraction of mRNA left", las = 1)
lines(tt, exp(-khat * tt), lwd = 3, col = oi["red"])
lines(tt, exp(-0.2 * tt), lwd = 2, lty = 2, col = oi["grey"])
lines(tt, exp(-1 * tt), lwd = 2, lty = 3, col = oi["grey"])
th <- log(2) / khat
segments(th, -0.1, th, 0.5, lty = 3); segments(-1, 0.5, th, 0.5, lty = 3)
text(th - 0.1, 0.45, sprintf("half-life %.2f h", th), adj = 1, cex = 1.3)
legend("topright", c(sprintf("best fit  k = %.3f / h", khat), "k = 0.2 / h", "k = 1 / h"),
       col = oi[c("red", "grey", "grey")], lwd = c(3, 2, 2), lty = c(1, 2, 3), cex = 1.3, bty = "n")
dev.off()

## ---- fig03: SSE(k) ------------------------------------------------------------

svg("fig03-decay-sse.svg", width = 9, height = 5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
kk <- seq(0, 3, length.out = 400)
plot(kk, ssev(kk), type = "n", ylim = c(0, 1.6),
     xlab = "decay rate k (per hour)", ylab = "SSE(k)", las = 1)
rect(kinfl, -1, 4, 2, col = "#F0E44255", border = NA)
lines(kk, ssev(kk), lwd = 3, col = oi["blue"])
points(khat, sse(khat), pch = 19, cex = 1.6, col = oi["red"])
text(khat + 0.08, 0.06, sprintf("minimum  k = %.4f", khat), adj = 0, cex = 1.3, col = oi["red"])
abline(v = kinfl, lty = 3, lwd = 2)
text(kinfl - 0.05, 1.45, "SSE'' > 0\n(convex)", adj = 1, cex = 1.3)
text(kinfl + 0.05, 1.45, "SSE'' < 0\n(concave)", adj = 0, cex = 1.3)
dev.off()

## ---- fig04: golden-section steps --------------------------------------------

svg("fig04-golden.svg", width = 10, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
kk <- seq(0, 1, length.out = 300)
plot(kk, ssev(kk), type = "l", lwd = 3, col = oi["blue"], ylim = c(-0.36, 0.4), xlim = c(-0.12, 1),
     xlab = "decay rate k (per hour)", ylab = "SSE(k)", axes = FALSE)
axis(1, at = seq(0, 1, 0.2)); axis(2, at = seq(0, 0.4, 0.1), las = 1); box()
a <- 0; b <- 1; cols <- oi[c("orange", "green", "purple", "red")]
for (i in 1:4) {
  x1 <- b - gr * (b - a); x2 <- a + gr * (b - a)
  yb <- -0.06 - 0.09 * (i - 1)
  segments(a, yb, b, yb, lwd = 4, col = cols[i])
  points(c(x1, x2), c(yb, yb), pch = 19, cex = 1.4, col = cols[i])
  text(-0.02, yb, paste0("step ", i), adj = 1, cex = 1.2)
  if (i < 4) {
    points(c(x1, x2), ssev(c(x1, x2)), pch = 19, cex = 1.3, col = cols[i])
    segments(c(x1, x2), ssev(c(x1, x2)), c(x1, x2), yb, lty = 3, col = cols[i])
  }
  if (sse(x1) < sse(x2)) b <- x2 else a <- x1
}
abline(h = 0, col = oi["grey"])
text(0.99, 0.37, "bars: bracket at each step\n(width x 0.618 per step)", adj = 1, cex = 1.2)
dev.off()

## ---- fig05: convergence -------------------------------------------------------

golden_trace <- function(f, a, b, nev) {
  xs <- numeric(0); fs <- numeric(0)
  ev <- function(x) { v <- f(x); xs <<- c(xs, x); fs <<- c(fs, v); v }
  x1 <- b - gr * (b - a); f1 <- ev(x1); x2 <- a + gr * (b - a); f2 <- ev(x2)
  while (length(xs) < nev) {
    if (f1 < f2) { b <- x2; x2 <- x1; f2 <- f1; x1 <- b - gr * (b - a); f1 <- ev(x1) }
    else { a <- x1; x1 <- x2; f1 <- f2; x2 <- a + gr * (b - a); f2 <- ev(x2) }
  }
  sapply(seq_along(xs), function(i) xs[which.min(fs[1:i])])
}
best_err <- function(xs, f) { fs <- sapply(xs, f); sapply(seq_along(xs), function(i) abs(xs[which.min(fs[1:i])] - khat)) }
ge <- abs(golden_trace(sse, 0, 1, 45) - khat)
bx <- numeric(0); invisible(optimize(function(k) { bx <<- c(bx, k); sse(k) }, c(0, 1), tol = 1e-12))
be <- best_err(bx, sse)
nk <- 0.2; for (i in 1:6) nk <- c(nk, nk[i] - d1(nk[i]) / d2(nk[i]))
ne <- abs(nk - khat)
floor_ <- 1e-16
svg("fig05-convergence.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(seq_along(ge), pmax(ge, floor_), log = "y", type = "o", pch = 19, lwd = 2, col = oi["orange"],
     ylim = c(1e-16, 1), xlim = c(1, 45), axes = FALSE,
     xlab = "function evaluations", ylab = "error in k")
axis(1); axis(2, at = 10^seq(-16, 0, 4), labels = c("1e-16", "1e-12", "1e-8", "1e-4", "1"), las = 1); box()
lines(seq_along(be), pmax(be, floor_), type = "o", pch = 17, lwd = 2, col = oi["blue"])
lines(2 * seq_along(ne) - 2 + 0.01, pmax(ne, floor_), type = "o", pch = 15, lwd = 2, col = oi["green"])
abline(h = sqrt(.Machine$double.eps) * khat, lty = 3, lwd = 2, col = oi["grey"])
text(45, 1.5e-8 * khat * 3, expression(sqrt(epsilon) %*% hat(k)), adj = 1, cex = 1.3)
legend("bottomleft", c("golden section (x 0.618 per evaluation)", "Brent (optimize, tol = 1e-12)",
                       "Newton on SSE' (2 evaluations per step: SSE', SSE'')"),
       col = oi[c("orange", "blue", "green")], pch = c(19, 17, 15), lwd = 2, cex = 1.15, bty = "n")
dev.off()

## ---- fig06: parabolic interpolation ---------------------------------------------

vertex <- function(a, b, c, fa, fb, fc)
  b - 0.5 * ((b - a)^2 * (fb - fc) - (b - c)^2 * (fb - fa)) / ((b - a) * (fb - fc) - (b - c) * (fb - fa))
parab <- function(xs, fs) { function(x) fs[1] * (x - xs[2]) * (x - xs[3]) / ((xs[1] - xs[2]) * (xs[1] - xs[3])) +
                                         fs[2] * (x - xs[1]) * (x - xs[3]) / ((xs[2] - xs[1]) * (xs[2] - xs[3])) +
                                         fs[3] * (x - xs[1]) * (x - xs[2]) / ((xs[3] - xs[1]) * (xs[3] - xs[2])) }
svg("fig06-parabola.svg", width = 12, height = 4.8)
par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
for (case in 1:2) {
  xs <- if (case == 1) c(0.2, 0.3, 0.5) else c(0, 0.5, 1)
  fs <- ssev(xs); v <- vertex(xs[1], xs[2], xs[3], fs[1], fs[2], fs[3]); P <- parab(xs, fs)
  kk <- if (case == 1) seq(0.15, 0.6, length.out = 300) else seq(-0.2, 1.1, length.out = 300)
  yl <- if (case == 1) c(-0.02, 0.15) else c(-0.6, 2.7)
  plot(kk, ssev(kk), type = "l", lwd = 3, col = oi["blue"], ylim = yl,
       xlab = "decay rate k", ylab = "SSE(k)",
       main = if (case == 1) "points close to the minimum" else "points far apart")
  lines(kk, P(kk), lwd = 2, lty = 2, col = oi["orange"])
  points(xs, fs, pch = 19, cex = 1.5, col = oi["orange"])
  abline(v = v, lty = 3, lwd = 2, col = oi["red"])
  text(v, yl[2] * 0.92, sprintf(" vertex %.3f", v), adj = 0, cex = 1.3, col = oi["red"])
  abline(v = khat, col = oi["grey"])
}
dev.off()

## ---- fig07: the flat bottom -------------------------------------------------

dd <- seq(-6e-8, 6e-8, length.out = 401)
svg("fig07-flat-bottom.svg", width = 12, height = 4.8)
par(mfrow = c(1, 2), mar = c(5, 6.5, 3, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
fv <- sapply(lhat + dd, nll) - nll(lhat)
plot(dd * 1e8, fv * 1e14, pch = 19, cex = 0.5, col = oi["blue"],
     xlab = expression(lambda - hat(lambda) ~ "(units of " * 10^-8 * ")"),
     ylab = expression(f(lambda) - f(hat(lambda)) ~ "(" * 10^-14 * ")"),
     main = "minus log-likelihood: a flat, noisy floor")
lines(dd * 1e8, 0.5 * 56.9186 * dd^2 * 1e14, lwd = 2, col = oi["red"])
plot(dd * 1e8, S(lhat + dd), pch = 19, cex = 0.5, col = oi["orange"],
     xlab = expression(lambda - hat(lambda) ~ "(units of " * 10^-8 * ")"),
     ylab = expression("score " * S(lambda)), main = "score: a clean sign change")
abline(h = 0, col = oi["grey"])
dev.off()

## ---- fig08: Newton's local parabola -------------------------------------------

svg("fig08-newton-min.svg", width = 12, height = 4.8)
par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
for (k0 in c(0.2, 1.0)) {
  kk <- if (k0 < 0.5) seq(0, 0.8, length.out = 300) else seq(0, 4.5, length.out = 300)
  q <- function(k) sse(k0) + d1(k0) * (k - k0) + 0.5 * d2(k0) * (k - k0)^2
  k1 <- k0 - d1(k0) / d2(k0)
  yl <- if (k0 < 0.5) c(0, 0.4) else c(0, 1.6)
  plot(kk, ssev(kk), type = "l", lwd = 3, col = oi["blue"], ylim = yl,
       xlab = "decay rate k", ylab = "SSE(k)",
       main = sprintf("start k0 = %.1f:  SSE'' %s 0", k0, if (d2(k0) > 0) ">" else "<"))
  lines(kk, q(kk), lwd = 2, lty = 2, col = oi["orange"])
  points(k0, sse(k0), pch = 19, cex = 1.6, col = oi["orange"])
  abline(v = k1, lty = 3, lwd = 2, col = oi["red"])
  text(k1, yl[2] * 0.93, sprintf(" k1 = %.2f", k1), adj = if (k1 > max(kk) * 0.7) 1.05 else 0, cex = 1.3, col = oi["red"])
}
dev.off()

## ---- fig09: circadian period ---------------------------------------------------

tc <- seq(0, 140, by = 4)
set.seed(15)
yc <- round(2 + 0.6 * cos(2 * pi * (tc - 6) / 24.5) + rnorm(36, sd = 0.3), 2)
fitT <- function(T) lm.fit(cbind(1, cos(2 * pi * tc / T), sin(2 * pi * tc / T)), yc)
rss <- function(T) sum(fitT(T)$residuals^2)
Tg <- seq(8, 48, by = 0.02); rg <- sapply(Tg, rss)
T1 <- optimize(rss, c(15, 35), tol = 1e-8)$minimum
T2 <- optimize(rss, c(10, 40), tol = 1e-8)$minimum
svg("fig09-circadian.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
plot(tc, yc, pch = 19, cex = 1.2, col = oi["grey"], xlab = "time (h)", ylab = "bioluminescence (a.u.)",
     main = "data and two fitted sine waves", ylim = c(0.8, 3.4))
ts <- seq(0, 140, length.out = 600)
for (j in 1:2) {
  T <- c(T1, T2)[j]; cf <- fitT(T)$coefficients
  lines(ts, cf[1] + cf[2] * cos(2 * pi * ts / T) + cf[3] * sin(2 * pi * ts / T), lwd = 2.5,
        col = oi[c("blue", "red")][j], lty = j)
}
legend("top", c(sprintf("T = %.2f h", T1), sprintf("T = %.2f h", T2)), col = oi[c("blue", "red")],
       lwd = 2.5, lty = 1:2, horiz = TRUE, cex = 1.2, bty = "n")
plot(Tg, rg, type = "l", lwd = 2.5, col = oi["blue"], xlab = "period T (h)", ylab = "RSS(T)",
     main = "residual sum of squares", ylim = c(0, 11.5), las = 1)
points(c(T1, T2), sapply(c(T1, T2), rss), pch = 19, cex = 1.6, col = oi[c("blue", "red")])
p0 <- c(10 + (1 - gr) * 30, 10 + gr * 30)
points(p0, sapply(p0, rss), pch = 4, cex = 1.8, lwd = 3, col = oi["orange"])
text(p0, sapply(p0, rss) - 0.7, sprintf("%.1f", p0), cex = 1.1, col = oi["orange"])
text(T1 + 1, 1.5, "global minimum", adj = 0, cex = 1.2, col = oi["blue"])
text(T2, 9.2, "local", cex = 1.2, col = oi["red"])
dev.off()
