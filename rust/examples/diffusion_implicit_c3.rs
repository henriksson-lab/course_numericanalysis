// C3: Crank-Nicolson oscillations after a sharp start, how to detect them,
// and two fixes (backward Euler start steps, smaller time step).

/// Thomas algorithm (a[0] and c[n-1] unused), as in Structure: symmetric, banded, sparse.
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

/// One backward Euler step; c includes the two (fixed) boundary values.
fn be_step(c: &mut [f64], r: f64) {
    let m = c.len() - 2;
    let mut d = c[1..=m].to_vec();
    d[0] += r * c[0];
    d[m - 1] += r * c[m + 1];
    let x = thomas(&vec![-r; m], &vec![1.0 + 2.0 * r; m], &vec![-r; m], &d);
    c[1..=m].copy_from_slice(&x);
}

/// One Crank-Nicolson step.
fn cn_step(c: &mut [f64], r: f64) {
    let m = c.len() - 2;
    let mut d: Vec<f64> = (1..=m).map(|i| c[i] + r / 2.0 * (c[i - 1] - 2.0 * c[i] + c[i + 1])).collect();
    d[0] += r / 2.0 * c[0];
    d[m - 1] += r / 2.0 * c[m + 1];
    let x = thomas(&vec![-r / 2.0; m], &vec![1.0 + r; m], &vec![-r / 2.0; m], &d);
    c[1..=m].copy_from_slice(&x);
}

/// Exact chamber profile (Fourier series of the previous unit).
fn exact(x: f64, t: f64, l: f64, dd: f64) -> f64 {
    let pi = std::f64::consts::PI;
    let mut s = 1.0 - x / l;
    for n in 1..=200 {
        let k = n as f64 * pi / l;
        s -= 2.0 / (n as f64 * pi) * (k * x).sin() * (-k * k * dd * t).exp();
    }
    s
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max)
}

fn main() {
    let (dd, l, n) = (1000.0, 500.0, 50usize);
    let dx = l / n as f64;
    let r = dd * 4.0 / (dx * dx);                      // dt = 4 s, r = 40
    let mut empty = vec![0.0; n + 1];
    empty[0] = 1.0;
    let mut c = empty.clone();
    let mut c1 = vec![0.0; 15];
    for k in 0..15 { cn_step(&mut c, r); c1[k] = c[1]; }
    println!("{:?}", c1.iter().map(|v| (v * 1000.0).round() / 1000.0).collect::<Vec<_>>());
    let inc: Vec<f64> = c1.windows(2).map(|w| w[1] - w[0]).collect();
    let flips = inc.windows(2).filter(|w| w[0].signum() != w[1].signum()).count();
    println!("{}", flips);                              // sign changes of the increments
    println!("{}", c.iter().cloned().fold(f64::MIN, f64::max) > 1.0);
    let mut cr = empty.clone();                         // fix 1: 4 BE steps of 1 s, then CN
    for _ in 0..4 { be_step(&mut cr, r / 4.0); }
    for _ in 0..14 { cn_step(&mut cr, r); }
    let mut cs = empty.clone();                         // fix 2: CN with dt = 0.4 s
    for _ in 0..150 { cn_step(&mut cs, r / 10.0); }
    let ce: Vec<f64> = (0..=n).map(|i| exact(i as f64 * dx, 60.0, l, dd)).collect();
    println!("plain {:.6e}  rannacher {:.6e}  small {:.6e}",
             max_abs_diff(&c, &ce), max_abs_diff(&cr, &ce), max_abs_diff(&cs, &ce));
}
