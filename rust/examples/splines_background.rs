// Background: a natural cubic spline through the noisy OD data, hand-written
// (second derivatives from a tridiagonal solve, then the cubic on each interval).

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

fn main() {
    let t: Vec<f64> = (0..41).map(|i| i as f64 * 0.25).collect();   // every 15 min
    let y = [0.040, 0.058, 0.076, 0.077, 0.108, 0.127, 0.152, 0.191, 0.200,
             0.263, 0.284, 0.327, 0.382, 0.446, 0.502, 0.556, 0.610, 0.672,
             0.749, 0.794, 0.839, 0.883, 0.934, 0.958, 1.004, 1.032, 1.076,
             1.097, 1.105, 1.111, 1.144, 1.155, 1.163, 1.171, 1.167, 1.182,
             1.193, 1.184, 1.177, 1.197, 1.199];
    let m = spline_m(&t, &y);
    println!("{}", spline_eval(&t, &y, &m, 3.1, 0));
    println!("{}", spline_eval(&t, &y, &m, 3.1, 1));               // its slope
}
