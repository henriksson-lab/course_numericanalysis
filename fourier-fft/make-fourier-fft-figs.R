# Figures for the fourier-fft/ unit ("Fourier series and the FFT").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd fourier-fft && Rscript make-fourier-fft-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). The timing figure (fig04) takes about a minute.
# Makes:
#   fig01-circadian-data.svg   a clock gene sampled every 2 h for 48 h
#   fig02-trig-interp.svg      trigonometric interpolant: periodic data vs a growth curve
#   fig03-fft-split.svg        the even/odd split of the FFT for n = 8
#   fig04-timing.svg           time of DFT matrix, hand-written FFT and fft() vs n
#   fig05-aliasing.svg         a 4.8-h rhythm sampled every 2 h and every 4 h
#   fig06-leakage.svg          spectrum of a 24-h and a 25.5-h rhythm, with/without Hann window
#   fig07-smoothing.svg        low-pass filtering a noisy luciferase trace, raw and detrended
#   fig08-blur.svg             image blur as convolution, and naive deconvolution

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## The running example: a clock-controlled gene, log2 expression, every 2 h for 48 h.
## Generated once with set.seed(1); 8 + 1.5 cos(2 pi (t - 6)/24) + N(0, 0.3^2), rounded.
x <- c(7.81, 8.81, 9.05, 9.98, 9.40, 8.50, 8.15, 7.47, 6.87, 6.41, 7.15, 7.37,
       7.81, 8.09, 9.64, 9.49, 9.29, 9.03, 8.25, 7.43, 6.98, 6.73, 6.72, 6.65)
t <- seq(0, 46, by = 2)

trig_interp <- function(x, T, tt) {
  n <- length(x); X <- fft(x) / n
  p <- Re(X[1]) + Re(X[n/2 + 1]) * cos(pi * n * tt / T)
  for (k in 1:(n/2 - 1))
    p <- p + 2 * Re(X[k + 1]) * cos(2 * pi * k * tt / T) -
             2 * Im(X[k + 1]) * sin(2 * pi * k * tt / T)
  p
}

## ---- fig01: the data -------------------------------------------------------
svg("fig01-circadian-data.svg", width = 11, height = 4.6)
par(mar = c(4.8, 5.2, 1, 1), cex.lab = 1.5, cex.axis = 1.35)
plot(t, x, pch = 19, col = oi["blue"], cex = 1.6, xaxt = "n",
     xlab = "time (h)", ylab = "log2 expression", ylim = c(6, 10.3))
axis(1, at = seq(0, 48, by = 6))
rect(c(12, 36), 5, c(24, 48), 11, col = "#00000012", border = NA)
lines(t, x, col = oi["blue"])
text(c(6, 18, 30, 42), 10.2, c("light", "dark", "light", "dark"), cex = 1.3,
     col = "grey30")
dev.off()

## ---- fig02: trig interpolation: periodic data vs non-periodic data ---------
tt <- seq(0, 47.9, by = 0.1)
p <- trig_interp(x, 48, tt)
X <- fft(x) / 24
fit2 <- Re(X[1]) + 2 * Re(X[3]) * cos(2 * pi * 2 * tt / 48) -
        2 * Im(X[3]) * sin(2 * pi * 2 * tt / 48)
## logistic growth curve (course running example) sampled every hour, 16 points
K <- 1.2; N0 <- 0.05; r <- 0.8
logis <- function(t) K / (1 + (K / N0 - 1) * exp(-r * t))
tg <- 0:15
og <- logis(tg)
ttg <- seq(0, 15.95, by = 0.05)
pg <- trig_interp(og, 16, ttg)

svg("fig02-trig-interp.svg", width = 13, height = 4.8)
par(mfrow = c(1, 2), mar = c(4.8, 5.2, 2.6, 1), cex.lab = 1.5, cex.axis = 1.35,
    cex.main = 1.5)
plot(tt, p, type = "l", lwd = 2, col = oi["blue"], xlab = "time (h)",
     ylab = "log2 expression", main = "clock gene: periodic", xaxt = "n",
     ylim = c(6, 10.5))
axis(1, at = seq(0, 48, by = 12))
lines(tt, fit2, lwd = 3, lty = 2, col = oi["orange"])
points(t, x, pch = 19, cex = 1.4)
legend("topright", c("interpolant (all 24 terms)", "keep k = 0 and 2 only"),
       col = oi[c("blue", "orange")], lwd = c(2, 3), lty = c(1, 2), cex = 1.15,
       bg = "white")
plot(ttg, pg, type = "l", lwd = 2, col = oi["red"], xlab = "time (h)",
     ylab = "OD", main = "growth curve: not periodic", ylim = c(-0.3, 1.5))
lines(ttg, logis(ttg), lwd = 2, lty = 2, col = oi["grey"])
points(tg, og, pch = 19, cex = 1.4)
legend("topleft", c("trig interpolant", "true curve"), col = oi[c("red", "grey")],
       lwd = 2, lty = c(1, 2), cex = 1.15, bg = "white")
dev.off()

## ---- fig03: even/odd split for n = 8 --------------------------------------
svg("fig03-fft-split.svg", width = 12, height = 5)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(NA, xlim = c(-2.4, 25), ylim = c(-0.3, 4.4), axes = FALSE, xlab = "", ylab = "")
rows <- list(list(0:7), list(c(0, 2, 4, 6), c(1, 3, 5, 7)),
             list(c(0, 4), c(2, 6), c(1, 5), c(3, 7)),
             as.list(c(0, 4, 2, 6, 1, 5, 3, 7)))
yy <- c(3.6, 2.4, 1.2, 0)
gcol <- c(oi["sky"], oi["yellow"], oi["green"], oi["purple"])
for (lev in 1:4) {
  pos <- 0
  for (g in rows[[lev]]) {
    for (v in g) {
      rect(pos * 2, yy[lev], pos * 2 + 1.6, yy[lev] + 0.7, col = gcol[lev], border = "white")
      text(pos * 2 + 0.8, yy[lev] + 0.35, bquote(x[.(v)]), cex = 1.5)
      pos <- pos + 1
    }
  }
  if (lev > 1) {
    gs <- length(rows[[lev]][[1]])
    for (b in seq(0, 7, by = gs))
      rect(b * 2 - 0.1, yy[lev] - 0.08, b * 2 + 2 * gs - 0.3, yy[lev] + 0.78,
           border = "grey30", lwd = 1.5)
  }
}
for (lev in 1:3) arrows(7.8, yy[lev] - 0.05, 7.8, yy[lev + 1] + 0.83, length = 0.12,
                        lwd = 2, col = "grey40")
text(-0.3, yy + 0.35, c("n = 8", "2 × 4", "4 × 2", "8 × 1"), adj = 1, cex = 1.4)
text(c(3.7, 11.7), 3.25, c("even samples -> E", "odd samples -> O"), cex = 1.2,
     col = "grey20")
arrows(17, 0.2, 17, 4.1, length = 0.15, lwd = 3, col = oi["red"])
text(17.6, 3.7, "combine upwards:", adj = 0, cex = 1.35, col = oi["red"])
text(17.6, 3.0, expression(X[k] == E[k] + w^k * O[k]), adj = 0, cex = 1.35, col = oi["red"])
text(17.6, 2.4, expression(X[k + n/2] == E[k] - w^k * O[k]), adj = 0, cex = 1.35,
     col = oi["red"])
text(17.6, 1.6, expression(w == e^{-2 * pi * i / n}), adj = 0, cex = 1.35, col = oi["red"])
text(17.6, 0.8, "n operations per level", adj = 0, cex = 1.35)
text(17.6, 0.25, expression(log[2] * 8 == 3 ~ "levels"), adj = 0, cex = 1.35)
dev.off()

## ---- fig04: timing ---------------------------------------------------------
my_fft <- function(x) {
  n <- length(x)
  if (n == 1) return(x)
  E <- my_fft(x[seq(1, n, by = 2)])
  O <- my_fft(x[seq(2, n, by = 2)])
  w <- exp(-2i * pi * (0:(n/2 - 1)) / n)
  c(E + w * O, E - w * O)
}
timeit <- function(f, x, reps) {
  system.time(for (i in 1:reps) f(x))[["elapsed"]] / reps
}
set.seed(1)
nm <- 2^(4:11); nr <- 2^(4:16); nf <- 2^(4:22)
tm <- sapply(nm, function(n) {
  xx <- rnorm(n); j <- 0:(n - 1)
  reps <- max(1, round(2^14 / n))
  system.time(for (i in 1:reps) {
    F <- exp(-2i * pi * outer(j, j) / n); F %*% xx })[["elapsed"]] / reps
})
tr <- sapply(nr, function(n) timeit(my_fft, rnorm(n), max(1, round(2^15 / n))))
tf <- sapply(nf, function(n) timeit(fft, rnorm(n), max(1, round(2^22 / n))))
tm <- pmax(tm, 1e-6); tr <- pmax(tr, 1e-6); tf <- pmax(tf, 1e-7)

svg("fig04-timing.svg", width = 9.5, height = 6)
par(mar = c(4.8, 6.4, 1, 1), cex.lab = 1.5, cex.axis = 1.3, mgp = c(4.2, 1, 0))
plot(nf, tf, log = "xy", type = "b", pch = 19, col = oi["green"], lwd = 2,
     xlim = range(nf), ylim = c(1e-7, 100), xaxt = "n", yaxt = "n",
     xlab = "", ylab = "time (seconds)")
mtext("n (signal length)", 1, line = 3, cex = 1.5)
axis(1, mgp = c(3, 1, 0), at = 10^(1:6), labels = c("10", "100", "1000", expression(10^4),
                                  expression(10^5), expression(10^6)))
axis(2, at = 10^seq(-6, 2, by = 2), labels = c(expression(10^-6), expression(10^-4),
     "0.01", "1", "100"), las = 1)
lines(nr, tr, type = "b", pch = 17, col = oi["orange"], lwd = 2)
lines(nm, tm, type = "b", pch = 15, col = oi["red"], lwd = 2)
nn <- 10^seq(3, 6, length.out = 20)
lines(nn, 3e-9 * nn^2, lty = 3, col = "grey40", lwd = 2)
text(2.5e4, 3e-9 * 1e8 * 2, expression("" %prop% n^2), cex = 1.4, col = "grey30")
legend("topleft", c("DFT matrix  F %*% x", "hand-written FFT (R)", "fft()"),
       col = oi[c("red", "orange", "green")], pch = c(15, 17, 19), lwd = 2,
       cex = 1.25, bg = "white")
dev.off()

## ---- fig05: aliasing -------------------------------------------------------
tc <- seq(0, 48, by = 0.05)
f48 <- function(t) cos(2 * pi * t / 4.8)
svg("fig05-aliasing.svg", width = 12, height = 4.6)
par(mar = c(4.8, 5.2, 1, 1), cex.lab = 1.5, cex.axis = 1.35)
plot(tc, f48(tc), type = "l", col = oi["grey"], lwd = 1.5, xaxt = "n",
     xlab = "time (h)", ylab = "expression", ylim = c(-1.25, 1.55))
axis(1, at = seq(0, 48, by = 4))
lines(tc, cos(2 * pi * tc / 24), col = oi["orange"], lwd = 3, lty = 2)
points(seq(0, 46, by = 2), f48(seq(0, 46, by = 2)), pch = 19, col = oi["blue"], cex = 1.5)
points(seq(0, 44, by = 4), f48(seq(0, 44, by = 4)), pch = 21, bg = oi["orange"],
       cex = 2.2, lwd = 1.5)
legend("top", horiz = TRUE, c("true rhythm, period 4.8 h", "every 2 h", "every 4 h",
       "24-h cosine through the 4-h samples"),
       col = c(oi["grey"], oi["blue"], "black", oi["orange"]),
       pt.bg = c(NA, NA, oi["orange"], NA), pch = c(NA, 19, 21, NA),
       lty = c(1, NA, NA, 2), lwd = c(1.5, NA, 1.5, 3), cex = 1.1, bty = "n")
dev.off()

## ---- fig06: leakage and the Hann window -----------------------------------
tl <- seq(0, 119.5, by = 0.5); nl <- length(tl)
hann <- 0.5 - 0.5 * cos(2 * pi * (0:(nl - 1)) / nl)
spec <- function(y, w = rep(1, length(y))) 2 * Mod(fft(w * y))[1:41] / sum(w)
s24 <- spec(cos(2 * pi * tl / 24))
s25 <- spec(cos(2 * pi * tl / 25.5))
s25h <- spec(cos(2 * pi * tl / 25.5), hann)
kk <- 0:40
svg("fig06-leakage.svg", width = 11, height = 5)
par(mar = c(4.8, 5.4, 3.2, 1), cex.lab = 1.5, cex.axis = 1.3)
plot(kk, pmax(s24, 1e-5), log = "y", type = "h", lwd = 6, col = oi["grey"],
     ylim = c(1e-5, 1.5), yaxt = "n", xlab = "k (cycles per 120 h)",
     ylab = "amplitude")
axis(2, at = 10^(-5:0), labels = c(expression(10^-5), expression(10^-4), "0.001",
     "0.01", "0.1", "1"), las = 1)
axis(3, at = c(1, 2, 3, 4, 5, 10, 20, 40), labels = c(120, 60, 40, 30, 24, 12, 6, 3),
     cex.axis = 1.1)
mtext("period (h)", 3, line = 2, cex = 1.3)
lines(kk + 0.12, s25, type = "b", pch = 19, col = oi["red"], lwd = 2)
lines(kk - 0.12, s25h, type = "b", pch = 17, col = oi["blue"], lwd = 2)
legend("topright", c("period 24 h (exactly 5 cycles)", "period 25.5 h, no window",
       "period 25.5 h, Hann window"), col = oi[c("grey", "red", "blue")],
       lwd = c(6, 2, 2), pch = c(NA, 19, 17), cex = 1.2, bg = "white")
dev.off()

## ---- fig07: smoothing a noisy luciferase trace -----------------------------
set.seed(3)
ts <- seq(0, 119.5, by = 0.5); ns <- length(ts)
truth <- 20 + 60 * exp(-ts / 40) + 30 * exp(-ts / 80) * cos(2 * pi * (ts - 4) / 24)
y <- truth + rnorm(ns, sd = 6)
lowpass <- function(y, kmax) {
  n <- length(y); Y <- fft(y)
  k <- c(0:(n/2), -(n/2 - 1):-1)
  Y[abs(k) > kmax] <- 0
  Re(fft(Y, inverse = TRUE)) / n
}
s1 <- lowpass(y, 15)
fit <- lm(y ~ ts)
s2 <- lowpass(resid(fit), 15) + fitted(fit)
svg("fig07-smoothing.svg", width = 13, height = 4.8)
par(mfrow = c(1, 2), mar = c(4.8, 5.2, 2.6, 1), cex.lab = 1.5, cex.axis = 1.35,
    cex.main = 1.45)
for (i in 1:2) {
  plot(ts, y, pch = 16, cex = 0.6, col = oi["grey"], xlab = "time (h)",
       ylab = "bioluminescence", xaxt = "n", ylim = c(0, 125),
       main = c("low-pass on the raw trace", "remove a straight line first")[i])
  axis(1, at = seq(0, 120, by = 24))
  lines(ts, truth, lwd = 2, lty = 2)
  lines(ts, list(s1, s2)[[i]], lwd = 3, col = oi[c("red", "blue")[i]])
}
legend("topright", c("noisy data", "true signal", "smoothed"), pch = c(16, NA, NA),
       lty = c(NA, 2, 1), lwd = c(NA, 2, 3), col = c(oi["grey"], "black", oi["blue"]),
       cex = 1.15, bg = "white")
dev.off()

## ---- fig08: image blur as convolution -------------------------------------
set.seed(5)
N <- 128
img <- matrix(0, N, N)
cx <- sample(10:118, 25); cy <- sample(10:118, 25)
for (i in 1:25) {
  d <- outer(1:N - cx[i], 1:N - cy[i], function(a, b) sqrt(a^2 + b^2))
  img[d <= 2.5] <- 1
}
d0 <- outer(c(0:(N/2), -(N/2 - 1):-1), c(0:(N/2), -(N/2 - 1):-1),
            function(a, b) a^2 + b^2)
psf <- exp(-d0 / (2 * 3^2)); psf <- psf / sum(psf)       # Gaussian, sd 3 px, centred at [1,1]
H <- fft(psf)
blur <- Re(fft(fft(img) * H, inverse = TRUE)) / N^2
noisy <- blur + matrix(rnorm(N^2, sd = 0.01), N)
deconv <- Re(fft(fft(noisy) / H, inverse = TRUE)) / N^2
shift <- function(m) m[c((N/2 + 1):N, 1:(N/2)), c((N/2 + 1):N, 1:(N/2))]
gr <- gray.colors(256, start = 0, end = 1)
show <- function(m, main) {
  image(1:nrow(m), 1:ncol(m), m, col = gr, axes = FALSE, xlab = "", ylab = "", asp = 1,
        main = main, useRaster = TRUE)
}
svg("fig08-blur.svg", width = 14, height = 4.2)
par(mfrow = c(1, 4), mar = c(0.5, 0.5, 3, 0.5), cex.main = 1.6)
show(img, "cells (true image)")
show(shift(psf)[44:84, 44:84] / max(psf), "point spread function")
show(blur, "blurred = image * PSF")
show(pmin(pmax(deconv, -1), 2), "divide by fft(PSF) after noise")
dev.off()
