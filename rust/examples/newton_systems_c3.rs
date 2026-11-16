// C3: Newton from 100 starts; collect the distinct steady states.
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
    let g: Vec<f64> = (0..10).map(|i| 0.2 + 0.4 * i as f64).collect(); // 10 values
    let mut roots: Vec<DVector<f64>> = Vec::new();
    for &v in &g {
        for &u in &g {
            let p = newton_sys(f, jac, &v2(u, v), 1e-12, 50).unwrap();
            roots.push(p.last().unwrap().clone());   // last iterate
        }
    }
    // equal to 6 decimals = same root; count how often each one is found
    let mut distinct: Vec<([f64; 2], usize)> = Vec::new();
    for r in &roots {
        let key = [(r[0] * 1e6).round() / 1e6, (r[1] * 1e6).round() / 1e6];
        match distinct.iter_mut().find(|(k, _)| *k == key) {
            Some((_, c)) => *c += 1,
            None => distinct.push((key, 1)),
        }
    }
    for (k, c) in &distinct {
        println!("{:?}  {:.2}", k, *c as f64 / roots.len() as f64);
    }
    let worst = roots.iter().map(|x| f(x).norm()).fold(0.0, f64::max);
    println!("{:e}", worst);                    // worst residual of all 100 runs
}
