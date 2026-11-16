// C4: Newton with a finite-difference Jacobian, for several step sizes h.
use nalgebra::{DMatrix, DVector};

const ALPHA: f64 = 3.0;
const N: f64 = 2.0;

fn f(x: &DVector<f64>) -> DVector<f64> {
    let (u, v) = (x[0], x[1]);
    DVector::from_vec(vec![ALPHA / (1.0 + v.powf(N)) - u,
                           ALPHA / (1.0 + u.powf(N)) - v])
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
    let xstar = v2((3.0 + 5f64.sqrt()) / 2.0, (3.0 - 5f64.sqrt()) / 2.0);
    for h in [f64::EPSILON.sqrt(), 1e-4, 1e-2, 0.3] {
        let jfd = |x: &DVector<f64>| fd_jacobian(f, x, h);
        let p = newton_sys(f, jfd, &v2(3.0, 0.5), 1e-12, 50).unwrap();
        let errs: Vec<String> = p[1..6].iter()
            .map(|x| format!("{:.0e}", (x - &xstar).norm())).collect();
        println!("h = {:.1e}: {} steps; errors {}", h, p.len() - 1, errs.join(" "));
    }
}
