// C4: method of lines on the chamber (N = 200). ode_solvers has no stiff solver:
// the explicit Dopri5 (= RK45) shows the cost of stiffness; the implicit side is our
// hand-written Crank-Nicolson with Thomas (4 backward Euler start steps).
use ode_solvers::{DVector, Dopri5, OutputType, System};

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

struct Chamber { dd: f64, dx: f64 }
impl System<f64, DVector<f64>> for Chamber {
    fn system(&self, _t: f64, c: &DVector<f64>, dc: &mut DVector<f64>) {
        let m = c.len();                                // interior points only
        let k = self.dd / (self.dx * self.dx);
        for i in 0..m {
            let left = if i == 0 { 1.0 } else { c[i - 1] };       // source 1
            let right = if i == m - 1 { 0.0 } else { c[i + 1] };  // sink 0
            dc[i] = k * (left - 2.0 * c[i] + right);
        }
    }
}

fn main() {
    let (dd, l, n) = (1000.0, 500.0, 200usize);
    let dx = l / n as f64;
    let ce: Vec<f64> = (0..=n).map(|i| exact(i as f64 * dx, 60.0, l, dd)).collect();

    // explicit Dopri5, stiffness test switched off
    let mut s = Dopri5::from_param(Chamber { dd, dx }, 0.0, 60.0, 60.0, DVector::zeros(n - 1),
                                   1e-6, 1e-6, 0.9, 0.04, 0.2, 10.0, 60.0, 0.0,
                                   1_000_000, u32::MAX, OutputType::Sparse);
    let st = s.integrate().unwrap();
    let y = s.y_out().last().unwrap();
    let err = (0..n - 1).map(|i| (y[i] - ce[i + 1]).abs()).fold(0.0, f64::max);
    println!("Dopri5       steps {:6}  f-evals {:6}  error {:.3e}", st.accepted_steps, st.num_eval, err);

    // hand-written Crank-Nicolson, dt = 0.25 s, after 4 backward Euler steps of 0.0625 s
    let mut c = vec![0.0; n + 1];
    c[0] = 1.0;
    let r = dd * 0.25 / (dx * dx);
    for _ in 0..4 { be_step(&mut c, r / 4.0); }
    for _ in 0..239 { cn_step(&mut c, r); }
    println!("CN by hand   steps    243  solves    243  error {:.3e}", max_abs_diff(&c, &ce));
}
