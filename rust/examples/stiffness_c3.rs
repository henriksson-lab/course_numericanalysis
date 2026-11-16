// C3: step counts and time on the enzyme model. Rust's ode_solvers has only explicit
// solvers (Dopri5, Dop853, Rk4); the stiff side is our hand-written backward Euler from C2.
// Run with --release for meaningful timings.
use nalgebra::{DMatrix, DVector};
use ode_solvers::{Dopri5, OutputType, System, Vector4};
use std::time::Instant;

const KF: f64 = 100.0; const KR: f64 = 100.0; const KCAT: f64 = 1.0;

struct Enzyme;
impl System<f64, Vector4<f64>> for Enzyme {
    fn system(&self, _t: f64, y: &Vector4<f64>, dy: &mut Vector4<f64>) {
        let b = KF * y[1] * y[0] - KR * y[2];
        dy[0] = -b; dy[1] = -b + KCAT * y[2]; dy[2] = b - KCAT * y[2]; dy[3] = KCAT * y[2];
    }
}

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
fn beuler(y0: &DVector<f64>, h: f64, n: usize) -> (DVector<f64>, usize) {  // (y, f evaluations)
    let (mut y, mut nf) = (y0.clone(), 0);
    let id = DMatrix::<f64>::identity(4, 4);
    for _ in 0..n {
        let yold = y.clone();
        for _ in 0..20 {
            let g = &y - &yold - h * f(&y); nf += 1;
            let delta = (&id - h * jac(&y)).lu().solve(&(-g)).unwrap();
            y += &delta;
            if delta.norm() < 1e-10 { break; }
        }
    }
    (y, nf)
}

fn main() {
    let y0 = Vector4::new(10.0, 1.0, 0.0, 0.0);
    // Dopri5 with the stiffness test switched off (n_stiff = u32::MAX)
    let t = Instant::now();
    let mut s = Dopri5::from_param(Enzyme, 0.0, 30.0, 30.0, y0, 1e-6, 1e-9, 0.9, 0.04,
                                   0.2, 10.0, 30.0, 0.0, 1_000_000, u32::MAX, OutputType::Sparse);
    let st = s.integrate().unwrap();
    println!("Dopri5  steps {:5}  f-evals {:6}  time {:.4} s  P(30) = {:.6}", st.accepted_steps,
             st.num_eval, t.elapsed().as_secs_f64(), s.y_out().last().unwrap()[3]);
    // hand-written backward Euler, fixed h = 0.01 s (no error control: order 1)
    let t = Instant::now();
    let (y, nf) = beuler(&DVector::from_vec(vec![10.0, 1.0, 0.0, 0.0]), 0.01, 3000);
    println!("beuler  steps  3000  f-evals {:6}  time {:.4} s  P(30) = {:.6}", nf,
             t.elapsed().as_secs_f64(), y[3]);
}
