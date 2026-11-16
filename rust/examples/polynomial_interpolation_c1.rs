// C1: three forms of the same interpolating polynomial.
use nalgebra::{DMatrix, DVector};

const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

fn lagrange(x: &[f64], y: &[f64], z: f64) -> f64 {
    let mut p = 0.0;
    for i in 0..x.len() {
        let mut li = 1.0;
        for j in 0..x.len() {
            if j != i {
                li *= (z - x[j]) / (x[i] - x[j]);
            }
        }
        p += li * y[i];
    }
    p
}

fn divdiff(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut a = y.to_vec();
    for k in 1..n {
        for i in (k..n).rev() {
            a[i] = (a[i] - a[i - 1]) / (x[i] - x[i - k]);
        }
    }
    a
}

fn newton_eval(x: &[f64], a: &[f64], z: f64) -> f64 {
    let n = a.len();
    let mut p = a[n - 1];
    for k in (0..n - 1).rev() {
        p = a[k] + (z - x[k]) * p;
    }
    p
}

fn main() {
    let t = [0.0, 2.0, 4.0, 6.0, 8.0];
    let y: Vec<f64> = t.iter().map(|&ti| od(ti)).collect();
    // (a) monomial coefficients from the Vandermonde system (nalgebra LU)
    let v = DMatrix::from_fn(5, 5, |i, j| t[i].powi(j as i32));
    let cf = v.lu().solve(&DVector::from_vec(y.clone())).unwrap();
    let a = divdiff(&t, &y);
    for z in [1.0f64, 3.0, 5.0, 7.0] {
        let p_mono: f64 = (0..5).map(|j| cf[j] * z.powi(j as i32)).sum();
        let p_lag = lagrange(&t, &y, z);
        let p_new = newton_eval(&t, &a, z);
        println!("{z:2.0}  {p_mono:.7}  {p_lag:.7}  {p_new:.7}  truth {:.7}  error {:+.4}",
                 od(z), p_lag - od(z));
    }
    // (b) check B1 and B2
    let (x3, y3) = ([0.0, 2.0, 4.0], [0.05, 0.21, 0.62]);
    println!("{:?}", [lagrange(&x3, &y3, 1.0), lagrange(&x3, &y3, 3.0)]);
    let (x4, y4) = ([0.0, 2.0, 4.0, 6.0], [0.05, 0.21, 0.62, 1.01]);
    let a4 = divdiff(&x4, &y4);
    println!("{:?}", a4);
    println!("{}", newton_eval(&x4, &a4, 3.0));
}
