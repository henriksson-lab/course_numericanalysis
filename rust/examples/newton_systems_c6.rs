// C6: failures (nearly singular J, NaN for a non-integer Hill coefficient)
// and a damped Newton that keeps concentrations non-negative.
use nalgebra::{DMatrix, DVector};

const ALPHA: f64 = 3.0;

// the Hill coefficient n is a parameter here, so we can change it
fn f(x: &DVector<f64>, n: f64) -> DVector<f64> {
    let (u, v) = (x[0], x[1]);
    DVector::from_vec(vec![ALPHA / (1.0 + v.powf(n)) - u,
                           ALPHA / (1.0 + u.powf(n)) - v])
}

fn jac(x: &DVector<f64>, n: f64) -> DMatrix<f64> {
    let (u, v) = (x[0], x[1]);
    DMatrix::from_row_slice(2, 2, &[
        -1.0, -ALPHA * n * v.powf(n - 1.0) / (1.0 + v.powf(n)).powi(2),
        -ALPHA * n * u.powf(n - 1.0) / (1.0 + u.powf(n)).powi(2), -1.0])
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

fn newton_damped(n: f64, x0: &DVector<f64>, tol: f64, maxit: usize) -> DVector<f64> {
    let mut x = x0.clone();
    for k in 1..=maxit {
        let delta = jac(&x, n).lu().solve(&(-f(&x, n))).expect("singular Jacobian");
        let mut t = 1.0;
        // halve the step until x stays >= 0 and ||F|| does not grow
        loop {
            let xt = &x + &delta * t;
            if xt.iter().all(|&c| c >= 0.0) && f(&xt, n).norm() <= f(&x, n).norm() {
                break;
            }
            t /= 2.0;
            if t < 1e-10 { panic!("damping failed"); }
        }
        x += &delta * t;
        println!("{} {} {:?}", k, t, x.as_slice());
        if (&delta * t).norm() < tol { break; }
    }
    x
}

fn main() {
    // (a) two nearby starts
    for s in [[1.9, 0.5], [2.0, 0.55]] {
        let x0 = v2(s[0], s[1]);
        let p = newton_sys(|x: &DVector<f64>| f(x, 2.0), |x: &DVector<f64>| jac(x, 2.0), &x0, 1e-12, 50).unwrap();
        let step1 = (&p[1] - &p[0]).norm();
        let minval = p.iter().flat_map(|x| x.iter().copied()).fold(f64::INFINITY, f64::min);
        let end = p.last().unwrap();
        println!("start {:?} | det J: {:.3} | first step length: {:.3} | min u, v on the way: {:.3} | end: {:.6} {:.6}",
                 s, jac(&x0, 2.0).determinant(), step1, minval, end[0], end[1]);
    }
    // (b) Hill coefficient 2.5: plain Newton
    let mut x = v2(2.0, 0.5);
    for _ in 0..2 {
        // with NaN entries the LU still "succeeds" and returns NaN
        x += jac(&x, 2.5).lu().solve(&(-f(&x, 2.5))).unwrap();
        println!("{:?}", x.as_slice());
    }
    // (c) damped Newton
    let xd = newton_damped(2.5, &v2(2.0, 0.5), 1e-12, 50);
    println!("{:?}", f(&xd, 2.5).as_slice());
}
