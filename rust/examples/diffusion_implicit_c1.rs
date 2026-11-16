// C1: backward Euler and Crank-Nicolson steps (hand-written), checked by hand (B1)
// and against the exact chamber solution at t = 60 s.

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
    let c0 = [1.0, 0.0, 0.0, 0.0, 0.0];             // B1: N = 4, r = 1
    let mut c1 = c0.to_vec();
    be_step(&mut c1, 1.0);
    println!("{:?}", c1.iter().map(|v| v * 21.0).collect::<Vec<_>>());
    let mut c2 = c1.clone();
    be_step(&mut c2, 1.0);
    println!("{:?}", c2);
    let mut c3 = c0.to_vec();
    cn_step(&mut c3, 1.0);
    println!("{:?}", c3);

    let (dd, l, n) = (1000.0, 500.0, 50usize);
    let dx = l / n as f64;
    let ce: Vec<f64> = (0..=n).map(|i| exact(i as f64 * dx, 60.0, l, dd)).collect();
    for dt in [4.0, 1.0, 0.25] {
        let r = dd * dt / (dx * dx);
        let mut cb = vec![0.0; n + 1];
        cb[0] = 1.0;
        let mut cc = cb.clone();
        for _ in 0..(60.0 / dt) as usize { be_step(&mut cb, r); cn_step(&mut cc, r); }
        println!("dt = {:4.2}  r = {:5.1}  BE error {:.2e}  CN error {:.2e}",
                 dt, r, max_abs_diff(&cb, &ce), max_abs_diff(&cc, &ce));
    }
}
