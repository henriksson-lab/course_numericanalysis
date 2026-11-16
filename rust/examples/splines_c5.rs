// C5: growth rate from the noisy OD reads: central differences, the slope of the
// interpolating spline, and the slope of a smoothing spline with lam chosen by GCV.
// Hand-written (nalgebra for the small dense solve); scipy's make_smoothing_spline
// solves the same problem.
use nalgebra::{DMatrix, DVector};

// Thomas algorithm (Structure: symmetric, banded, sparse); a[0] and c[n-1] unused
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec());
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    x
}

// second derivatives M of the natural spline (M[0] = M[n-1] = 0)
fn spline_m(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let s: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / h[i]).collect();
    let a: Vec<f64> = h[..n - 2].to_vec();                       // h_{i-1}
    let c: Vec<f64> = h[1..].to_vec();                           // h_i
    let b: Vec<f64> = (0..n - 2).map(|i| 2.0 * (a[i] + c[i])).collect();
    let d: Vec<f64> = (0..n - 2).map(|i| 6.0 * (s[i + 1] - s[i])).collect();
    let mut m = vec![0.0];
    m.extend(thomas(&a, &b, &c, &d));
    m.push(0.0);
    m
}

// value (deriv = 0) or slope (deriv = 1) of the spline at t
fn spline_eval(x: &[f64], y: &[f64], m: &[f64], t: f64, deriv: u8) -> f64 {
    let i = x.partition_point(|&xi| xi <= t).clamp(1, x.len() - 1) - 1;
    let h = x[i + 1] - x[i];
    let (a, b) = ((x[i + 1] - t) / h, (t - x[i]) / h);
    if deriv == 0 {
        a * y[i] + b * y[i + 1] + ((a * a * a - a) * m[i] + (b * b * b - b) * m[i + 1]) * h * h / 6.0
    } else {
        (y[i + 1] - y[i]) / h + (-(3.0 * a * a - 1.0) * m[i] + (3.0 * b * b - 1.0) * m[i + 1]) * h / 6.0
    }
}


// Smoothing spline (Reinsch form): minimise sum (y_i - g_i)^2 + lam * int S''^2.
// With Q (n x n-2) and R (n-2 x n-2, tridiagonal) built from the spacings,
// solve (R + lam Q^T Q) gamma = Q^T y; fitted values g = y - lam Q gamma;
// gamma are the second derivatives at the interior knots (natural spline).
struct Fit { g: Vec<f64>, m: Vec<f64>, df: f64 }

fn smooth_fit(x: &[f64], y: &[f64], lam: f64) -> Fit {
    let n = x.len();
    let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let mut q = DMatrix::<f64>::zeros(n, n - 2);
    let mut r = DMatrix::<f64>::zeros(n - 2, n - 2);
    for j in 0..n - 2 {                            // column j <-> interior knot j+1
        q[(j, j)] = 1.0 / h[j];
        q[(j + 1, j)] = -1.0 / h[j] - 1.0 / h[j + 1];
        q[(j + 2, j)] = 1.0 / h[j + 1];
        r[(j, j)] = (h[j] + h[j + 1]) / 3.0;
        if j + 1 < n - 2 { r[(j, j + 1)] = h[j + 1] / 6.0; r[(j + 1, j)] = h[j + 1] / 6.0; }
    }
    let yv = DVector::from_column_slice(y);
    let binv = (&r + lam * q.transpose() * &q).try_inverse().unwrap();
    let gamma = &binv * (q.transpose() * &yv);
    let g = &yv - lam * &q * &gamma;
    let s = lam * &q * &binv * q.transpose();      // hat matrix A = I - S
    let df: f64 = (0..n).map(|i| 1.0 - s[(i, i)]).sum();  // trace of A
    let mut m = vec![0.0];
    m.extend(gamma.iter());
    m.push(0.0);
    Fit { g: g.iter().cloned().collect(), m, df }
}

fn gcv(x: &[f64], y: &[f64], lam: f64) -> f64 {
    let f = smooth_fit(x, y, lam);
    let n = x.len() as f64;
    let rss: f64 = y.iter().zip(&f.g).map(|(a, b)| (a - b).powi(2)).sum();
    n * rss / (n - f.df).powi(2)
}

// minimise a function of log10(lam): coarse grid, then golden section
fn minimise_log(fun: &dyn Fn(f64) -> f64, lo: f64, hi: f64) -> f64 {
    let grid: Vec<f64> = (0..=90).map(|i| lo + (hi - lo) * i as f64 / 90.0).collect();
    let k = (0..grid.len()).min_by(|&a, &b| fun(grid[a]).partial_cmp(&fun(grid[b])).unwrap()).unwrap();
    let (mut a, mut b) = (grid[k.saturating_sub(1)], grid[(k + 1).min(90)]);
    let phi = (5f64.sqrt() - 1.0) / 2.0;
    while b - a > 1e-8 {
        let (c, d) = (b - phi * (b - a), a + phi * (b - a));
        if fun(c) < fun(d) { b = d } else { a = c }
    }
    0.5 * (a + b)
}

fn main() {
    let t: Vec<f64> = (0..41).map(|i| i as f64 * 0.25).collect();
    let y = [0.040, 0.058, 0.076, 0.077, 0.108, 0.127, 0.152, 0.191, 0.200,
             0.263, 0.284, 0.327, 0.382, 0.446, 0.502, 0.556, 0.610, 0.672,
             0.749, 0.794, 0.839, 0.883, 0.934, 0.958, 1.004, 1.032, 1.076,
             1.097, 1.105, 1.111, 1.144, 1.155, 1.163, 1.171, 1.167, 1.182,
             1.193, 1.184, 1.177, 1.197, 1.199];
    let (k, n0, r) = (1.2, 0.05, 0.8);
    let od = |t: f64| k / (1.0 + (k / n0 - 1.0) * (-r * t).exp());
    let rate = |t: f64| r * od(t) * (1.0 - od(t) / k);
    let rms = |e: &[f64]| (e.iter().map(|v| v * v).sum::<f64>() / e.len() as f64).sqrt();
    let n = t.len();
    let cd: Vec<f64> = (1..n - 1).map(|i| (y[i + 1] - y[i - 1]) / 0.5 - rate(t[i])).collect();
    println!("{}", rms(&cd));                                  // central differences
    let m = spline_m(&t, &y);                                   // interpolating spline
    let e: Vec<f64> = t.iter().map(|&ti| spline_eval(&t, &y, &m, ti, 1) - rate(ti)).collect();
    println!("{}", rms(&e));
    let l10 = minimise_log(&|l| gcv(&t, &y, 10f64.powf(l)), -6.0, 3.0);
    let fit = smooth_fit(&t, &y, 10f64.powf(l10));             // smoothing spline, GCV
    println!("lam {:.4}  df {:.4}", 10f64.powf(l10), fit.df);
    let e: Vec<f64> = t.iter().map(|&ti| spline_eval(&t, &fit.g, &fit.m, ti, 1) - rate(ti)).collect();
    println!("{}", rms(&e));
    let (mut best, mut at) = (f64::NEG_INFINITY, 0.0);
    for i in 0..=1000 {
        let ti = i as f64 / 100.0;
        let d = spline_eval(&t, &fit.g, &fit.m, ti, 1);
        if d > best { best = d; at = ti; }
    }
    println!("{:.7} {:.2}", best, at);                          // true: 0.24 at 3.92 h
    for lam in [10.0, 0.005] {
        let f = smooth_fit(&t, &y, lam);
        let e: Vec<f64> = t.iter().map(|&ti| spline_eval(&t, &f.g, &f.m, ti, 1) - rate(ti)).collect();
        println!("lam {:5}  df {:.2}  rms {}", lam, f.df, rms(&e));
    }
}
