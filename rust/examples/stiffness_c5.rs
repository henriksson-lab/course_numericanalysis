// C5: order of backward Euler and BDF2 on the enzyme model, from t = 0.1 s to 10.1 s.
// Reference values from Dopri5 at tight tolerance (stiffness test off).
use nalgebra::{DMatrix, DVector};
use ode_solvers::{Dopri5, OutputType, System, Vector4};

const KF: f64 = 100.0; const KR: f64 = 100.0; const KCAT: f64 = 1.0;

fn f(y: &DVector<f64>) -> DVector<f64> {
    let b = KF * y[1] * y[0] - KR * y[2];
    DVector::from_vec(vec![-b, -b + KCAT * y[2], b - KCAT * y[2], KCAT * y[2]])
}
fn jac(y: &DVector<f64>) -> DMatrix<f64> {
    DMatrix::from_row_slice(4, 4, &[
        -KF * y[1], -KF * y[0],  KR,         0.0,
        -KF * y[1], -KF * y[0],  KR + KCAT,  0.0,
         KF * y[1],  KF * y[0], -KR - KCAT,  0.0,
         0.0,        0.0,        KCAT,       0.0])
}

// one implicit step: solve y - c*h*f(y) = rhs by Newton, starting from y_start
fn implicit_solve(rhs: &DVector<f64>, c: f64, h: f64, y_start: &DVector<f64>) -> DVector<f64> {
    let id = DMatrix::<f64>::identity(4, 4);
    let mut y = y_start.clone();
    for _ in 0..20 {
        let g = &y - c * h * f(&y) - rhs;
        let delta = (&id - c * h * jac(&y)).lu().solve(&(-g)).unwrap();
        y += &delta;
        if delta.norm() < 1e-10 { break; }
    }
    y
}
fn beuler(y0: &DVector<f64>, h: f64, n: usize) -> DVector<f64> {
    let mut y = y0.clone();
    for _ in 0..n { y = implicit_solve(&y, 1.0, h, &y); }          // y - h f(y) = y_n
    y
}
fn bdf2(y0: &DVector<f64>, y1: &DVector<f64>, h: f64, n: usize) -> DVector<f64> {
    let (mut yprev, mut y) = (y0.clone(), y1.clone());
    for _ in 1..n {                                                 // y - 2/3 h f(y) = 4/3 y_n - 1/3 y_{n-1}
        let ynew = implicit_solve(&(4.0 / 3.0 * &y - 1.0 / 3.0 * &yprev), 2.0 / 3.0, h, &y);
        yprev = y; y = ynew;
    }
    y
}

struct Enzyme;
impl System<f64, Vector4<f64>> for Enzyme {
    fn system(&self, _t: f64, y: &Vector4<f64>, dy: &mut Vector4<f64>) {
        let b = KF * y[1] * y[0] - KR * y[2];
        dy[0] = -b; dy[1] = -b + KCAT * y[2]; dy[2] = b - KCAT * y[2]; dy[3] = KCAT * y[2];
    }
}
fn reference(t: f64) -> DVector<f64> {
    let mut s = Dopri5::from_param(Enzyme, 0.0, t, t, Vector4::new(10.0, 1.0, 0.0, 0.0),
                                   1e-12, 1e-14, 0.9, 0.04, 0.2, 10.0, t, 0.0, 10_000_000,
                                   u32::MAX, OutputType::Sparse);
    s.integrate().unwrap();
    DVector::from_column_slice(s.y_out().last().unwrap().as_slice())
}

fn main() {
    let ys = reference(0.1);                                        // after the binding transient
    let pref = reference(10.1)[3];
    for h in [0.2, 0.1, 0.05, 0.02, 0.01] {
        let n = (10.0 / h).round() as usize;
        let e1 = (beuler(&ys, h, n)[3] - pref).abs();
        let e2 = (bdf2(&ys, &reference(0.1 + h), h, n)[3] - pref).abs();
        println!("h = {:4.2}   backward Euler {:.3e}   BDF2 {:.3e}", h, e1, e2);
    }
}
