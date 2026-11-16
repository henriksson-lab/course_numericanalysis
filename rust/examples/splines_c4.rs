// C4: overshoot of a cubic spline through a gene time course, and two fixes:
// a monotone (PCHIP, as scipy's PchipInterpolator) interpolant and a spline of log2.
// All hand-written; the natural spline stands in for the library default.

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

// PCHIP slopes (Fritsch-Butland weighted harmonic mean, scipy's end rule)
fn pchip_slopes(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let h: Vec<f64> = (0..n - 1).map(|i| x[i + 1] - x[i]).collect();
    let s: Vec<f64> = (0..n - 1).map(|i| (y[i + 1] - y[i]) / h[i]).collect();
    let mut d = vec![0.0; n];
    for k in 1..n - 1 {
        if s[k - 1] * s[k] > 0.0 {                 // same sign: harmonic mean
            let (w1, w2) = (2.0 * h[k] + h[k - 1], h[k] + 2.0 * h[k - 1]);
            d[k] = (w1 + w2) / (w1 / s[k - 1] + w2 / s[k]);
        }                                          // local extremum: slope 0
    }
    let end = |h0: f64, h1: f64, s0: f64, s1: f64| {
        let d0 = ((2.0 * h0 + h1) * s0 - h0 * s1) / (h0 + h1);
        if d0.signum() != s0.signum() { 0.0 }
        else if s0.signum() != s1.signum() && d0.abs() > 3.0 * s0.abs() { 3.0 * s0 }
        else { d0 }
    };
    d[0] = end(h[0], h[1], s[0], s[1]);
    d[n - 1] = end(h[n - 2], h[n - 3], s[n - 2], s[n - 3]);
    d
}

// cubic Hermite interpolant with values y and slopes d at the knots
fn hermite_eval(x: &[f64], y: &[f64], d: &[f64], t: f64) -> f64 {
    let i = x.partition_point(|&xi| xi <= t).clamp(1, x.len() - 1) - 1;
    let h = x[i + 1] - x[i];
    let u = (t - x[i]) / h;
    let (h00, h10) = ((1.0 + 2.0 * u) * (1.0 - u).powi(2), u * (1.0 - u).powi(2));
    let (h01, h11) = (u * u * (3.0 - 2.0 * u), u * u * (u - 1.0));
    h00 * y[i] + h10 * h * d[i] + h01 * y[i + 1] + h11 * h * d[i + 1]
}

fn main() {
    let tg: [f64; 8] = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0];   // hours after stimulation
    let eg: [f64; 8] = [1.0, 40.0, 32.0, 6.0, 1.0, 0.5, 0.5, 0.5];   // TPM
    let lg: Vec<f64> = eg.iter().map(|v| v.log2()).collect();
    let (m, ml, dp) = (spline_m(&tg, &eg), spline_m(&tg, &lg), pchip_slopes(&tg, &eg));
    let fits: [(&str, Box<dyn Fn(f64) -> f64>); 3] = [
        ("natural", Box::new(|t| spline_eval(&tg, &eg, &m, t, 0))),
        ("monotone", Box::new(|t| hermite_eval(&tg, &eg, &dp, t))),
        ("log2", Box::new(|t| 2f64.powf(spline_eval(&tg, &lg, &ml, t, 0)))),
    ];
    let tt: Vec<f64> = (0..=2400).map(|i| i as f64 / 100.0).collect();
    for (name, f) in fits.iter() {
        let v: Vec<f64> = tt.iter().map(|&t| f(t)).collect();
        let (imin, vmin) = v.iter().enumerate().fold((0, f64::INFINITY),
            |acc, (i, &vi)| if vi < acc.1 - 1e-12 { (i, vi) } else { acc });
        let vmax = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        println!("{:9} at_3h {:7.3}  min {:7.3}  t_min {:5.2}  max {:6.3}",
                 name, f(3.0), vmin, tt[imin], vmax);
    }
}

// Output:
// natural   at_3h  -1.083  min  -1.100  t_min  3.06  max 41.725
// monotone  at_3h   2.311  min   0.500  t_min  8.00  max 40.000
// log2      at_3h   2.034  min   0.496  t_min 17.07  max 52.756
