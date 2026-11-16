// C2: stencil weights from a small linear system (nalgebra).
use nalgebra::{DMatrix, DVector};
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}
fn factorial(k: usize) -> f64 {
    (1..=k).map(|i| i as f64).product()
}

fn fd_weights(s: &[f64], m: usize) -> Vec<f64> {
    let n = s.len();
    let a = DMatrix::from_fn(n, n, |k, j| s[j].powi(k as i32) / factorial(k));
    let mut b = DVector::zeros(n);
    b[m] = 1.0; // pick out the m-th derivative
    a.lu().solve(&b).unwrap().iter().cloned().collect()
}

fn main() {
    println!("{:?}", fd_weights(&[-1.0, 0.0, 1.0], 1));
    println!("{:?}", fd_weights(&[-1.0, 0.0, 1.0], 2));
    println!("{:?}", fd_weights(&[0.0, 1.0, 2.0], 1));
    let s5 = [-2.0, -1.0, 0.0, 1.0, 2.0];
    let w5 = fd_weights(&s5, 1);
    let w12: Vec<f64> = w5.iter().map(|w| (w * 12.0 * 1e10).round() / 1e10 + 0.0).collect();
    println!("{:?}", w12); // five-point, times 12
    let exact = R * od(3.0) * (1.0 - od(3.0) / K);
    let d5 = |h: f64| s5.iter().zip(&w5).map(|(s, w)| w * od(3.0 + s * h)).sum::<f64>() / h;
    let (e1, e2) = (d5(0.1) - exact, d5(0.05) - exact);
    println!("{:e} {:e} {}", e1, e2, e1 / e2);
}
