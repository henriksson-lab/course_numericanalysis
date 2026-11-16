// C6: smoothing a gene along pseudotime; lam by GCV and by leave-one-out CV,
// a too-rough and a too-stiff fit, and the lag-1 residual correlation check.
// Hand-written smoothing spline (nalgebra for the dense solve), as in C5.
use nalgebra::{DMatrix, DVector};

// scientific notation as R and Python print it: 1.114e-04 rather than Rust's 1.114e-4
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
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
struct Fit { g: Vec<f64>, m: Vec<f64>, hat_diag: Vec<f64>, df: f64 }

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
    let hat_diag: Vec<f64> = (0..n).map(|i| 1.0 - s[(i, i)]).collect();
    let df = hat_diag.iter().sum();
    let mut m = vec![0.0];
    m.extend(gamma.iter());
    m.push(0.0);
    Fit { g: g.iter().cloned().collect(), m, hat_diag, df }
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

// leave-one-out CV: for a linear smoother the LOO residual is r_i / (1 - A_ii)
fn loocv(x: &[f64], y: &[f64], lam: f64) -> f64 {
    let f = smooth_fit(x, y, lam);
    let n = x.len() as f64;
    (0..x.len()).map(|i| ((y[i] - f.g[i]) / (1.0 - f.hat_diag[i])).powi(2)).sum::<f64>() / n
}

fn main() {
    let pt = [0.013, 0.062, 0.108, 0.126, 0.177, 0.186, 0.202, 0.206, 0.212, 0.266,
              0.267, 0.340, 0.372, 0.380, 0.382, 0.384, 0.386, 0.411, 0.482, 0.494,
              0.498, 0.573, 0.600, 0.629, 0.652, 0.661, 0.668, 0.687, 0.718, 0.724,
              0.770, 0.777, 0.794, 0.827, 0.870, 0.898, 0.908, 0.935, 0.945, 0.992];
    let ex = [1.32, 1.27, 1.03, 0.30, 1.22, 0.98, 0.95, 0.49, 0.84, 1.16,
              1.49, 1.06, 1.33, 1.21, 0.76, 1.11, 1.12, 1.39, 2.46, 2.48,
              2.20, 2.85, 3.24, 3.12, 2.53, 2.45, 2.76, 2.70, 2.04, 2.32,
              1.69, 1.28, 1.49, 0.81, 1.58, 1.73, 0.90, 0.65, 1.21, 0.95];
    let lam_gcv = 10f64.powf(minimise_log(&|l| gcv(&pt, &ex, 10f64.powf(l)), -7.0, -1.0));
    let lam_cv = 10f64.powf(minimise_log(&|l| loocv(&pt, &ex, 10f64.powf(l)), -7.0, -1.0));
    for (name, lam) in [("gcv", lam_gcv), ("loocv", lam_cv), ("rough", 1e-6), ("stiff", 1e-2)] {
        let f = smooth_fit(&pt, &ex, lam);
        let (mut height, mut at) = (f64::NEG_INFINITY, 0.0);
        for i in 0..=1000 {
            let g = i as f64 / 1000.0;
            let v = spline_eval(&pt, &f.g, &f.m, g, 0);
            if v > height { height = v; at = g; }
        }
        let e: Vec<f64> = ex.iter().zip(&f.g).map(|(a, b)| a - b).collect();
        let (a, b) = (&e[1..], &e[..e.len() - 1]);
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let (ma, mb) = (mean(a), mean(b));
        let cov: f64 = a.iter().zip(b).map(|(p, q)| (p - ma) * (q - mb)).sum();
        let va: f64 = a.iter().map(|p| (p - ma).powi(2)).sum();
        let vb: f64 = b.iter().map(|q| (q - mb).powi(2)).sum();
        println!("{:6} lam {}  df {:.4}  at {:.3}  height {:.4}  lag1 {:.4}",
                 name, sci(lam, 3), f.df, at, height, cov / (va * vb).sqrt());
    }
}
