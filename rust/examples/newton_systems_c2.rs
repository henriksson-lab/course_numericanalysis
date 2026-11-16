// C2: Newton for the toggle switch; quadratic convergence.
use nalgebra::{DMatrix, DVector};

const ALPHA: f64 = 3.0;
const N: f64 = 2.0;

fn f(x: &DVector<f64>) -> DVector<f64> {
    let (u, v) = (x[0], x[1]);
    DVector::from_vec(vec![ALPHA / (1.0 + v.powf(N)) - u,
                           ALPHA / (1.0 + u.powf(N)) - v])
}

fn jac(x: &DVector<f64>) -> DMatrix<f64> {
    let (u, v) = (x[0], x[1]);
    DMatrix::from_row_slice(2, 2, &[
        -1.0, -ALPHA * N * v.powf(N - 1.0) / (1.0 + v.powf(N)).powi(2),
        -ALPHA * N * u.powf(N - 1.0) / (1.0 + u.powf(N)).powi(2), -1.0])
}

// returns the path x_0, x_1, ...; None if J is exactly singular at some step
fn newton_sys(f: impl Fn(&DVector<f64>) -> DVector<f64>,
              jac: impl Fn(&DVector<f64>) -> DMatrix<f64>,
              x0: &DVector<f64>, tol: f64, maxit: usize) -> Option<Vec<DVector<f64>>> {
    let mut x = x0.clone();
    let mut path = vec![x.clone()];
    for _ in 0..maxit {
        let delta = jac(&x).lu().solve(&(-f(&x)))?; // one linear solve per step
        x += &delta;
        path.push(x.clone());
        if delta.norm() < tol {                      // norm() = sqrt(sum of squares)
            return Some(path);
        }
    }
    eprintln!("warning: no convergence");
    Some(path)
}

fn v2(a: f64, b: f64) -> DVector<f64> { DVector::from_vec(vec![a, b]) }

fn main() {
    let xstar = v2((3.0 + 5f64.sqrt()) / 2.0, (3.0 - 5f64.sqrt()) / 2.0); // exact, from B5
    let p = newton_sys(f, jac, &v2(3.0, 0.0), 1e-12, 50).unwrap();
    for x in &p[0..3] {
        println!("{:.6} {:.6}", x[0], x[1]);   // compare with B4
    }
    let p = newton_sys(f, jac, &v2(3.0, 0.5), 1e-12, 50).unwrap();
    let err: Vec<f64> = p.iter().map(|x| (x - &xstar).norm()).collect();
    for (k, x) in p.iter().enumerate() {
        println!("{}  {:.2e}  {:.2e}", k, err[k], f(x).norm());
    }
    let ratio: Vec<f64> = (0..4).map(|k| err[k + 1] / err[k].powi(2)).collect();
    println!("{:?}", ratio);                    // quadratic: roughly constant
}
