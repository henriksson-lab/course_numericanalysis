# Figures for the floating-point/ unit ("How computers store numbers").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd floating-point && Rscript make-floating-point-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Makes:
#   fig01-double-bits.svg    the 64 bits of an IEEE double
#   fig02-toy-floats.svg     a toy float system on the number line
#   fig03-spacing.svg        gap between neighbouring doubles vs |x|
#   fig04-range.svg          normal / subnormal / underflow / overflow ranges
#   fig05-underflow.svg      likelihood of an alignment: product vs log-sum

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- fig01: the 64 bits of a double ---------------------------------------
## Bits of 0.1 laid out as sign | exponent (11) | fraction (52).

bits_of <- function(x) {
  # IEEE bits of a double, most significant first (works on little-endian)
  raw <- rev(writeBin(x, raw(), endian = "little"))
  as.integer(unlist(lapply(raw, function(b) rev(as.integer(rawToBits(b))))))
}
b <- bits_of(0.1)

svg("fig01-double-bits.svg", width = 12, height = 3.2)
par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(NA, xlim = c(0, 64), ylim = c(-1.6, 2.0), axes = FALSE, xlab = "", ylab = "")
col <- c(oi["red"], rep(oi["sky"], 11), rep(oi["yellow"], 52))
rect(0:63, 0, 1:64, 1, col = col, border = "white", lwd = 1.5)
text(0:63 + 0.5, 0.5, b, cex = 1.0)
text(0, -0.35, "sign: 1 bit", cex = 1.3, adj = 0)
text(6.5, 1.45, "exponent: 11 bits", cex = 1.3)
text(38, 1.45, "fraction: 52 bits", cex = 1.3)
text(32, -1.1, expression(0.1 %~~% +1.6 %*% 2^-4 ~~~~~~
     "value" == (-1)^s %*% "1.f" %*% 2^{e - 1023}), cex = 1.5)
dev.off()

## ---- fig02: a toy float system on the number line --------------------------
## Numbers 1.f x 2^e with 3 fraction bits and e in -1..2.

toy <- sort(unique(as.vector(outer(1 + (0:7) / 8, 2^(-1:2)))))
svg("fig02-toy-floats.svg", width = 12, height = 3.6)
par(mar = c(4.5, 1, 1, 1))
plot(NA, xlim = c(0, 8), ylim = c(-0.6, 1.3), axes = FALSE,
     xlab = "", ylab = "")
abline(h = 0, col = oi["grey"])
ecol <- c(oi["blue"], oi["green"], oi["orange"], oi["purple"])
for (k in 0:3) {
  v <- (1 + (0:7) / 8) * 2^(k - 1)
  segments(v, -0.25, v, 0.25, col = ecol[k + 1], lwd = 3)
  text(mean(range(v)), c(0.95, 0.5)[k %% 2 + 1],
       paste0("[", 2^(k - 1), ", ", 2^k, "): gap ", c("1/16", "1/8", "1/4", "1/2")[k + 1]),
       col = ecol[k + 1], cex = 1.3)
}
segments(0, -0.25, 0, 0.25, lwd = 3)
axis(1, at = 0:8, cex.axis = 1.4, pos = -0.4)
mtext("3 fraction bits: 8 numbers per power of 2, gap doubles each time",
      side = 1, line = 3, cex = 1.3)
dev.off()

## ---- fig03: spacing between neighbouring doubles ---------------------------

eps <- .Machine$double.eps
x <- 10^seq(-3, 17, length.out = 4000)
gap <- 2^(floor(log2(x)) - 52)
svg("fig03-spacing.svg", width = 9, height = 5.5)
par(mar = c(5, 7.5, 1, 1), cex.lab = 1.4, cex.axis = 1.3, mgp = c(3.2, 1, 0))
plot(x, gap, log = "xy", type = "s", lwd = 3, col = oi["blue"],
     xlab = "|x|", ylab = "", axes = FALSE)
lines(x, eps * x, lty = 2, lwd = 2, col = oi["red"])
axis(1, at = 10^seq(-3, 17, 4), labels = parse(text = paste0("10^", seq(-3, 17, 4))))
axis(2, at = 10^seq(-19, 1, 4), labels = parse(text = paste0("10^", seq(-19, 1, 4))), las = 1)
title(ylab = "gap to the next double", line = 5.5)
box()
points(c(1, 2^53), c(eps, 2), pch = 19, cex = 1.6, col = oi["orange"])
text(2, eps / 50, expression("x = 1: gap = " * epsilon %~~% 2.2 %*% 10^-16),
     adj = 0, cex = 1.2)
text(2^53, 2, expression(2^53 * ": gap = 2"), pos = 2, cex = 1.2)
legend("topleft", c("actual gap", expression(epsilon %*% "|x|")),
       lty = c(1, 2), lwd = c(3, 2), col = c(oi["blue"], oi["red"]),
       bty = "n", cex = 1.3)
dev.off()

## ---- fig04: the range of doubles -------------------------------------------

svg("fig04-range.svg", width = 12, height = 3.4)
par(mar = c(4.5, 1, 1, 1))
plot(NA, xlim = c(-340, 340), ylim = c(0, 1.3), axes = FALSE,
     xlab = "", ylab = "")
rect(-340, 0.2, -323.3, 0.8, col = oi["grey"], border = NA)
rect(-323.3, 0.2, -307.7, 0.8, col = oi["yellow"], border = NA)
rect(-307.7, 0.2, 308.25, 0.8, col = oi["sky"], border = NA)
rect(308.25, 0.2, 340, 0.8, col = oi["red"], border = NA)
text(-331.5, 0.5, "0", cex = 1.3)
text(-326, 1.05, "subnormal (down to 4.9e-324)", cex = 1.3, adj = 0)
segments(-315, 0.82, -315, 0.95)
text(0, 0.5, "normal doubles: 15-16 significant digits", cex = 1.4)
text(324, 0.5, "Inf", cex = 1.3)
text(340, 1.05, "overflow above 1.8e308", cex = 1.3, adj = 1)
axis(1, at = c(-308, -200, -100, 0, 100, 200, 308),
     cex.axis = 1.2, pos = 0.15)
mtext(expression(log[10] ~ "|x|"), side = 1, line = 3, cex = 1.4)
dev.off()

## ---- fig05: likelihood of an alignment underflows --------------------------
## Same simulation as on the slides: 10 000 per-site likelihoods.

set.seed(1)
p <- runif(10000, min = 0.01, max = 0.5)     # per-site likelihoods
naive <- log10(cumprod(p))                   # -Inf once the product is 0
logsum <- cumsum(log10(p))
n <- seq_along(p)
first0 <- which(cumprod(p) == 0)[1]

svg("fig05-underflow.svg", width = 9, height = 5.5)
par(mar = c(5, 6, 1, 1), cex.lab = 1.4, cex.axis = 1.3)
plot(n[1:1000], logsum[1:1000], type = "l", lwd = 3, col = oi["blue"],
     xlab = "number of sites multiplied", ylab = expression(log[10] ~ "likelihood"),
     ylim = c(-700, 0), las = 1)
ok <- is.finite(naive[1:1000])
lines(n[1:1000][ok], naive[1:1000][ok], lwd = 6, col = oi["orange"], lty = 1)
lines(n[1:1000], logsum[1:1000], lwd = 2, col = oi["blue"])
abline(h = log10(.Machine$double.xmin), lty = 2, col = oi["grey"], lwd = 2)
abline(h = log10(2^-1074), lty = 3, col = oi["grey"], lwd = 2)
points(first0, -700, pch = 4, cex = 2.5, lwd = 4, col = oi["red"])
text(first0, -700, paste0("site ", first0, ": prod(p) = 0"), pos = 4,
     cex = 1.3, col = oi["red"])
text(1000, -280, "smallest normal double", pos = 2, cex = 1.1)
text(1000, -345, "smallest subnormal", pos = 2, cex = 1.1)
legend("topright", c("prod(p) in doubles", "sum(log10(p))"),
       lwd = c(6, 2), col = c(oi["orange"], oi["blue"]), bty = "n", cex = 1.3)
dev.off()
