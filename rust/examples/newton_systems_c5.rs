// C5: no library root finder in our crate list, so this is the hand-written
// Newton with a finite-difference Jacobian (what multiroot does by default),
// counting the evaluations of F.
use nalgebra::{DMatrix, DVector};
use std::cell::Cell;

const ALPHA: f64 = 3.0;
const N: f64 = 2.0;

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

fn fd_jacobian(f: impl Fn(&DVector<f64>) -> DVector<f64>, x: &DVector<f64>, h: f64) -> DMatrix<f64> {
    let fx = f(x);
    let mut jm = DMatrix::zeros(fx.len(), x.len());
    for j in 0..x.len() {
        let hj = h * x[j].abs().max(1.0);
        let mut xh = x.clone();
        xh[j] += hj;
        jm.set_column(j, &((f(&xh) - &fx) / hj));
    }
    jm
}

fn main() {
    let nfev = Cell::new(0);                    // a counter the closure may change
    let f = |x: &DVector<f64>| {
        nfev.set(nfev.get() + 1);
        let (u, v) = (x[0], x[1]);
        DVector::from_vec(vec![ALPHA / (1.0 + v.powf(N)) - u,
                               ALPHA / (1.0 + u.powf(N)) - v])
    };
    let jfd = |x: &DVector<f64>| fd_jacobian(&f, x, f64::EPSILON.sqrt());
    let xstar = v2((3.0 + 5f64.sqrt()) / 2.0, (3.0 - 5f64.sqrt()) / 2.0);
    let p = newton_sys(&f, jfd, &v2(3.0, 0.5), 1e-12, 50).unwrap();
    let x = p.last().unwrap();
    println!("root {:?}, steps {}, F evaluations {}", x.as_slice(), p.len() - 1, nfev.get());
    println!("residual {:?}", f(x).as_slice());
    println!("error {:?}", (x - &xstar).as_slice());
    for s in [[0.5, 3.0], [1.0, 1.0], [1.0, 0.2]] {
        let p = newton_sys(&f, |x: &DVector<f64>| fd_jacobian(&f, x, 1.5e-8), &v2(s[0], s[1]), 1e-12, 50).unwrap();
        println!("{:?}", p.last().unwrap().as_slice());
    }
}
