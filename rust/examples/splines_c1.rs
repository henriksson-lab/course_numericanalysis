// C1: max interpolation error vs spacing h on the exact logistic OD curve,
// piecewise linear and natural cubic spline (both hand-written).

// scientific notation as R and Python print it: 3.04e-02 rather than Rust's 3.04e-2
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

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

// piecewise linear interpolation (R approx, numpy interp)
fn lin_eval(x: &[f64], y: &[f64], t: f64) -> f64 {
    let i = x.partition_point(|&xi| xi <= t).clamp(1, x.len() - 1) - 1;
    y[i] + (y[i + 1] - y[i]) * (t - x[i]) / (x[i + 1] - x[i])
}

fn main() {
    let (k, n0, r) = (1.2, 0.05, 0.8);
    let od = |t: f64| k / (1.0 + (k / n0 - 1.0) * (-r * t).exp());
    let tt: Vec<f64> = (0..=1000).map(|i| i as f64 / 100.0).collect();   // fine grid
    for h in [2.0, 1.0, 0.5, 0.25, 0.125] {
        let m_pts = (10.0 / h) as usize;
        let x: Vec<f64> = (0..=m_pts).map(|i| i as f64 * h).collect();
        let y: Vec<f64> = x.iter().map(|&xi| od(xi)).collect();
        let m = spline_m(&x, &y);
        let (mut e_lin, mut e_nat) = (0.0_f64, 0.0_f64);
        for &t in &tt {
            e_lin = e_lin.max((lin_eval(&x, &y, t) - od(t)).abs());
            e_nat = e_nat.max((spline_eval(&x, &y, &m, t, 0) - od(t)).abs());
        }
        println!("h = {:5.3}  linear {}  natural {}", h, sci(e_lin, 2), sci(e_nat, 2));
    }
}

// Output:
// h = 2.000  linear 3.04e-02  natural 7.03e-03
// h = 1.000  linear 9.10e-03  natural 1.31e-03
// h = 0.500  linear 2.30e-03  natural 3.41e-04
// h = 0.250  linear 5.75e-04  natural 8.58e-05
// h = 0.125  linear 1.44e-04  natural 2.15e-05
