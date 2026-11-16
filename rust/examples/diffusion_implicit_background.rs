// Background: Thomas algorithm, one backward Euler step and one Crank-Nicolson
// step for the chamber (hand-written; no crate needed).

/// Thomas algorithm: a = below the diagonal (a[0] unused), b = diagonal,
/// c = above the diagonal (c[n-1] unused), d = right-hand side.
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

/// One backward Euler step; c includes the two boundary values (kept fixed).
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

fn main() {
    let (dd, l, n) = (1000.0, 500.0, 50);
    let dx = l / n as f64;
    let mut c = vec![0.0; n + 1];
    c[0] = 1.0;                                // empty channel, source 1, sink 0
    let r = dd * 4.0 / (dx * dx);              // dt = 4 s
    for _ in 0..15 { be_step(&mut c, r); }     // to t = 60 s
    let show = |c: &[f64]| {                   // x = 0, 100, ..., 500
        let s: Vec<String> = (0..=n).step_by(10).map(|i| format!("{:.4}", c[i])).collect();
        println!("[{}]", s.join(", "));
    };
    show(&c);
    let mut c2 = vec![0.0; n + 1];
    c2[0] = 1.0;
    for _ in 0..15 { cn_step(&mut c2, r); }
    show(&c2);
}
