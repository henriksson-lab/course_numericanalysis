# Figures for the newton/ unit ("Newton and secant methods").
#
# Re-run by hand from this folder (it is NOT run by `quarto render`):
#   cd newton && Rscript make-newton-figs.R
#
# Base R graphics only (grDevices::svg). Colours are Okabe-Ito
# (colour-blind safe). Makes:
#   fig01-two-site.svg       two-site dose-response curve and its c50
#   fig02-newton-secant.svg  Newton tangents and secant chords, side by side
#   fig03-convergence.svg    error vs iteration: bisection, secant, Newton
#   fig04-score.svg          score of the zero-truncated Poisson; bad Newton start
#   fig05-flat.svg           Newton from a high dose jumps to a negative dose
#   fig06-cycle.svg          Newton cycles on a Michaelis-Menten curve in log dose

oi <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73",
        yellow = "#F0E442", blue = "#0072B2", red = "#D55E00",
        purple = "#CC79A7", grey = "#999999")

## ---- the running examples (same code as on the slides) ---------------------

R  <- function(c) 0.3 * c / (0.1 + c) + 0.7 * c^2 / (25 + c^2)  # two-site response
f  <- function(c) R(c) - 0.5
df <- function(c) 0.03 / (0.1 + c)^2 + 35 * c / (25 + c^2)^2
c50 <- 3.261305265055212

x  <- rep(1:5, c(52, 30, 12, 5, 1))     # UMI counts in cells where the gene is seen
n  <- length(x); sx <- sum(x)
S  <- function(l) sx / l - n / (-expm1(-l))                     # score
dS <- function(l) -sx / l^2 + n * exp(-l) / expm1(-l)^2
lhat <- 1.218469525830130

newton_path <- function(f, df, x0, k) {
  xs <- x0
  for (i in 1:k) xs <- c(xs, xs[length(xs)] - f(xs[length(xs)]) / df(xs[length(xs)]))
  xs
}

## ---- fig01: two-site dose-response ------------------------------------------

svg("fig01-two-site.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
cc <- 10^seq(-3, 2, length.out = 400)
plot(cc, R(cc), type = "l", log = "x", lwd = 3, col = oi["blue"], ylim = c(0, 1),
     xlab = expression("dose " * c ~ "(" * mu * "M)"), ylab = "response R(c)", axes = FALSE)
axis(1, at = 10^(-3:2), labels = c("0.001", "0.01", "0.1", "1", "10", "100"))
axis(2, las = 1); box()
lines(cc, 0.3 * cc / (0.1 + cc), lty = 2, lwd = 2, col = oi["orange"])
lines(cc, 0.7 * cc^2 / (25 + cc^2), lty = 2, lwd = 2, col = oi["green"])
abline(h = 0.5, col = oi["grey"], lwd = 1.5)
segments(c50, -0.1, c50, 0.5, lty = 3, lwd = 2)
points(c50, 0.5, pch = 19, cex = 1.6, col = oi["red"])
text(c50 * 1.15, 0.43, expression(c[50] %~~% 3.26), adj = 0, cex = 1.4, col = oi["red"])
legend("topleft", c("total response", "high-affinity site (30 %)", "low-affinity site (70 %)"),
       col = oi[c("blue", "orange", "green")], lty = c(1, 2, 2), lwd = c(3, 2, 2),
       cex = 1.25, bty = "n")
dev.off()

## ---- fig02: Newton tangents vs secant chords --------------------------------

svg("fig02-newton-secant.svg", width = 12, height = 5)
par(mfrow = c(1, 2), mar = c(5, 5, 3, 1), cex.lab = 1.5, cex.axis = 1.3, cex.main = 1.5)
cc <- seq(0.5, 8, length.out = 300)
lab <- function(x, y, i, col) text(x, y, bquote(c[.(i)]), cex = 1.4, col = col)
## Newton from c0 = 7
plot(cc, f(cc), type = "l", lwd = 3, col = oi["blue"], xlab = expression("dose " * c),
     ylab = "f(c) = R(c) - 0.5", main = "Newton: follow the tangent", ylim = c(-0.25, 0.3))
abline(h = 0, col = oi["grey"])
p <- newton_path(f, df, 7, 2)
for (i in 1:2) {
  segments(p[i], f(p[i]), p[i + 1], 0, col = oi["red"], lwd = 2)
  segments(p[i + 1], 0, p[i + 1], f(p[i + 1]), lty = 3, col = oi["red"])
}
points(p, f(p), pch = 19, col = oi["red"], cex = 1.3)
for (i in 0:2) lab(p[i + 1], f(p[i + 1]) + c(0.03, -0.03, 0.03)[i + 1], i, oi["red"])
## secant from c0 = 7, c1 = 6
plot(cc, f(cc), type = "l", lwd = 3, col = oi["blue"], xlab = expression("dose " * c),
     ylab = "f(c) = R(c) - 0.5", main = "Secant: follow the chord", ylim = c(-0.25, 0.3))
abline(h = 0, col = oi["grey"])
s <- c(7, 6); for (i in 1:2) { k <- length(s); s <- c(s, s[k] - f(s[k]) * (s[k] - s[k - 1]) / (f(s[k]) - f(s[k - 1]))) }
for (i in 1:2) {
  lines(c(s[i], s[i + 2]), c(f(s[i]), 0), col = oi["green"], lwd = 2)
  segments(s[i + 2], 0, s[i + 2], f(s[i + 2]), lty = 3, col = oi["green"])
}
points(s, f(s), pch = 19, col = oi["green"], cex = 1.3)
for (i in 0:3) lab(s[i + 1], f(s[i + 1]) + c(0.03, 0.03, -0.03, 0.03)[i + 1], i, oi["green"])
dev.off()

## ---- fig03: convergence ------------------------------------------------------

lo <- 1; hi <- 10; eb <- c()
for (k in 1:40) { m <- (lo + hi) / 2; eb <- c(eb, abs(m - c50)); if (f(lo) * f(m) <= 0) hi <- m else lo <- m }
s <- c(1, 2); for (k in 1:8) { j <- length(s); if (f(s[j]) == f(s[j - 1])) break
  s <- c(s, s[j] - f(s[j]) * (s[j] - s[j - 1]) / (f(s[j]) - f(s[j - 1]))) }
es <- abs(s[-(1:2)] - c50)
en <- abs(newton_path(f, df, 1, 6)[-1] - c50)
fl <- 1e-17                              # plot exact zeros at the floor
svg("fig03-convergence.svg", width = 9, height = 5.5)
par(mar = c(5, 6.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3, mgp = c(3.8, 1, 0))
plot(seq_along(eb), pmax(eb, fl), log = "y", type = "o", pch = 15, lwd = 2, col = oi["grey"],
     xlim = c(1, 40), ylim = c(1e-17, 10), xlab = "iteration k", ylab = expression("|" * c[k] - c[50] * "|"),
     yaxt = "n")
axis(2, at = 10^seq(-16, 0, 4), labels = parse(text = paste0("10^", seq(-16, 0, 4))), las = 1)
lines(seq_along(es), pmax(es, fl), type = "o", pch = 17, lwd = 2, col = oi["green"])
lines(seq_along(en), pmax(en, fl), type = "o", pch = 19, lwd = 2, col = oi["red"])
abline(h = .Machine$double.eps * c50, lty = 3)
text(12, 3e-15, expression(epsilon %.% c[50]), adj = 0, cex = 1.3)
legend("bottomright", c("bisection on [1, 10]   (order 1)", "secant from 1, 2   (order 1.6)",
                       "Newton from 1   (order 2)"),
       col = oi[c("grey", "green", "red")], pch = c(15, 17, 19), lwd = 2, cex = 1.25, bg = "white")
dev.off()

## ---- fig04: score of the zero-truncated Poisson -----------------------------

svg("fig04-score.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
ll <- seq(0.4, 4, length.out = 300)
plot(ll, S(ll), type = "l", lwd = 3, col = oi["blue"], ylim = c(-80, 120),
     xlim = c(-0.8, 4), xlab = expression("rate " * lambda), ylab = expression("score " * S(lambda)))
abline(h = 0, col = oi["grey"]); abline(v = 0, col = oi["grey"], lty = 3)
rect(-2, -200, 0, 300, col = adjustcolor(oi["grey"], 0.2), border = NA)
text(-0.4, 100, expression(lambda < 0), cex = 1.3)
points(lhat, 0, pch = 19, col = oi["green"], cex = 1.6)
text(lhat + 0.05, 12, expression(hat(lambda) == 1.218), adj = 0, cex = 1.3, col = oi["green"])
l0 <- 3; l1 <- l0 - S(l0) / dS(l0)
segments(l0, S(l0), l1, 0, col = oi["red"], lwd = 2.5)
points(c(l0, l1), c(S(l0), 0), pch = 19, col = oi["red"], cex = 1.4)
text(l0, S(l0) - 10, expression(lambda[0] == 3), cex = 1.3, col = oi["red"])
text(l1, -12, expression(lambda[1] == -0.47), cex = 1.3, col = oi["red"])
abline(v = mean(x), lty = 2, col = oi["orange"], lwd = 2)
text(mean(x) + 0.05, 80, expression(bar(x) == 1.73), adj = 0, cex = 1.3, col = oi["orange"])
dev.off()

## ---- fig05: Newton from a high dose jumps to a negative dose -----------------

svg("fig05-flat.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
cc <- seq(-8, 12, length.out = 600)
cc <- cc[abs(cc + 0.1) > 0.15]           # skip the pole of c/(0.1 + c) at c = -0.1
fc <- f(cc); fc[abs(cc + 0.1) < 0.3] <- NA
plot(cc, fc, type = "l", lwd = 3, col = oi["blue"], ylim = c(-0.6, 0.6),
     xlab = expression("dose " * c ~ "(" * mu * "M)"), ylab = "f(c) = R(c) - 0.5")
rect(-20, -2, 0, 2, col = adjustcolor(oi["grey"], 0.2), border = NA)
text(-4, 0.5, "negative dose", cex = 1.3)
abline(h = 0, col = oi["grey"])
p <- newton_path(f, df, 10, 6)
segments(p[1], f(p[1]), p[2], 0, col = oi["red"], lwd = 2.5)
points(p[1], f(p[1]), pch = 19, col = oi["red"], cex = 1.4)
points(p[2], 0, pch = 19, col = oi["red"], cex = 1.4)
text(p[1], f(p[1]) + 0.07, expression(c[0] == 10), cex = 1.3, col = oi["red"])
text(p[2], -0.07, expression(c[1] == -5.73), cex = 1.3, col = oi["red"])
points(p[7], 0, pch = 4, lwd = 3, cex = 2, col = oi["red"])
segments(p[7], 0, -6, -0.3, col = oi["red"], lty = 3)
text(-7.8, -0.36, "Newton converges\nto c = -3.05", cex = 1.2, col = oi["red"], adj = 0)
points(c50, 0, pch = 19, cex = 1.6, col = oi["green"])
text(c50 + 0.6, -0.06, expression(c[50] == 3.26), adj = 0, cex = 1.3, col = oi["green"])
dev.off()

## ---- fig06: Newton cycles on Michaelis-Menten in log dose --------------------

us <- 2.177318984965307                  # root of sinh(u) = 2u
svg("fig06-cycle.svg", width = 9, height = 5.5)
par(mar = c(5, 5.5, 1, 1), cex.lab = 1.5, cex.axis = 1.3)
uu <- seq(-4, 4, length.out = 300)
g  <- function(u) 1 / (1 + exp(-u)) - 0.5
plot(uu, g(uu), type = "l", lwd = 3, col = oi["blue"], ylim = c(-0.55, 0.55),
     xlab = expression(u == ln(S / K[m])), ylab = expression(v / V[max] - 0.5))
abline(h = 0, v = 0, col = oi["grey"])
segments(us, g(us), -us, 0, col = oi["red"], lwd = 2.5)
segments(-us, g(-us), us, 0, col = oi["red"], lwd = 2.5)
segments(c(-us, us), 0, c(-us, us), g(c(-us, us)), lty = 3, col = oi["red"], lwd = 2)
points(c(-us, us), g(c(-us, us)), pch = 19, col = oi["red"], cex = 1.4)
text(us, g(us) + 0.07, expression(S[0] == 8.82 * K[m]), cex = 1.3, col = oi["red"])
text(-us, g(-us) - 0.07, expression(S[1] == K[m] / 8.82), cex = 1.3, col = oi["red"])
text(0.15, -0.45, expression("root: " * S == K[m]), adj = 0, cex = 1.3, col = oi["green"])
points(0, 0, pch = 19, col = oi["green"], cex = 1.5)
dev.off()
