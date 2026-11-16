// C3: natural vs clamped end conditions on 11 exact OD values (h = 1 h), hand-written.
// Rust has no not-a-knot / fmm library spline here, so the "default" row is missing.

// scientific notation as R and Python print it: 1.3134e-03 rather than Rust's 1.3134e-3
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

// clamped spline: S'(x_0) = d0, S'(x_{n-1}) = dn; full n x n tridiagonal system
fn clamped_m(x: &[f64], y: &[f64], d0: f64, dn: f64) -> Vec<f64> {
    let n = x.len();
    let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let s: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / h[i]).collect();
    let mut a = vec![0.0]; a.extend(&h);                  // a[i] = h_{i-1}
    let mut c = h.clone(); c.push(0.0);                   // c[i] = h_i
    let mut b = vec![2.0 * h[0]];
    b.extend((1..n - 1).map(|i| 2.0 * (h[i - 1] + h[i])));
    b.push(2.0 * h[n - 2]);
    let mut d = vec![6.0 * (s[0] - d0)];
    d.extend((1..n - 1).map(|i| 6.0 * (s[i] - s[i - 1])));
    d.push(6.0 * (dn - s[n - 2]));
    thomas(&a, &b, &c, &d)
}

fn main() {
    let (k, n0, r) = (1.2, 0.05, 0.8);
    let od = |t: f64| k / (1.0 + (k / n0 - 1.0) * (-r * t).exp());
    let rate = |t: f64| r * od(t) * (1.0 - od(t) / k);
    let x: Vec<f64> = (0..=10).map(|i| i as f64).collect();
    let y: Vec<f64> = x.iter().map(|&xi| od(xi)).collect();
    let tt: Vec<f64> = (0..=1000).map(|i| i as f64 / 100.0).collect();
    let cases = [("natural", spline_m(&x, &y)),
                 ("clamped", clamped_m(&x, &y, rate(0.0), rate(10.0))),
                 ("flat   ", clamped_m(&x, &y, 0.0, 0.0))];
    for (name, m) in cases.iter() {
        let (mut all, mut early) = (0.0_f64, 0.0_f64);
        for &t in &tt {
            let e = (spline_eval(&x, &y, m, t, 0) - od(t)).abs();
            all = all.max(e);
            if t <= 1.0 { early = early.max(e); }
        }
        println!("{} {} {}", name, sci(all, 4), sci(early, 4));
    }
}
