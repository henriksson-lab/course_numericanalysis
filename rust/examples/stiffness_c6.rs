// C6: Robertson's problem. (a) Dopri5 to t = 40: default settings stop with "stiffness
// detected"; with the test off it finishes, slowly. (b) Hand-written backward Euler with
// Newton on a log-spaced time grid all the way to t = 1e11. Run with --release for timings.
use nalgebra::{DMatrix, DVector};
use ode_solvers::{Dopri5, OutputType, System, Vector3};
use std::time::Instant;

struct Rober;
impl System<f64, Vector3<f64>> for Rober {
    fn system(&self, _t: f64, y: &Vector3<f64>, dy: &mut Vector3<f64>) {
        let (r1, r2, r3) = (0.04 * y[0], 3e7 * y[1] * y[1], 1e4 * y[1] * y[2]);
        dy[0] = -r1 + r3; dy[1] = r1 - r2 - r3; dy[2] = r2;
    }
}
fn f(y: &DVector<f64>) -> DVector<f64> {
    let (r1, r2, r3) = (0.04 * y[0], 3e7 * y[1] * y[1], 1e4 * y[1] * y[2]);
    DVector::from_vec(vec![-r1 + r3, r1 - r2 - r3, r2])
}
fn jac(y: &DVector<f64>) -> DMatrix<f64> {
    DMatrix::from_row_slice(3, 3, &[
        -0.04, 1e4 * y[2],                 1e4 * y[1],
         0.04, -6e7 * y[1] - 1e4 * y[2], -1e4 * y[1],
         0.0,   6e7 * y[1],                0.0])
}

fn main() {
    let y0 = Vector3::new(1.0, 0.0, 0.0);
    // (a) explicit Dopri5
    let mut s = Dopri5::new(Rober, 0.0, 40.0, 0.0, y0, 1e-6, 1e-10);
    if let Err(e) = s.integrate() { println!("Dopri5 (default): {}", e); }
    let t = Instant::now();
    let mut s = Dopri5::from_param(Rober, 0.0, 40.0, 0.0, y0, 1e-6, 1e-10, 0.9, 0.04, 0.2,
                                   10.0, 40.0, 0.0, 10_000_000, u32::MAX, OutputType::Sparse);
    let st = s.integrate().unwrap();
    println!("Dopri5 (test off): steps {}  f-evals {}  time {:.3} s  y2(40) = {:.6e}",
             st.accepted_steps, st.num_eval, t.elapsed().as_secs_f64(), s.y_out().last().unwrap()[1]);

    // (b) backward Euler on t = 0, 1e-6, ..., 1e11 with m points per decade
    for m in [10, 100] {
        let t = Instant::now();
        let times: Vec<f64> = std::iter::once(0.0)
            .chain((0..=17 * m).map(|k| 10f64.powf(-6.0 + k as f64 / m as f64))).collect();
        let mut y = DVector::from_vec(vec![1.0, 0.0, 0.0]);
        let id = DMatrix::<f64>::identity(3, 3);
        let mut out = Vec::new();
        for k in 1..times.len() {
            let h = times[k] - times[k - 1];
            let yold = y.clone();
            for _ in 0..20 {
                let g = &y - &yold - h * f(&y);
                let delta = (&id - h * jac(&y)).lu().solve(&(-g)).unwrap();
                y += &delta;
                if delta.norm() < 1e-12 * (1.0 + y.norm()) { break; }
            }
            if (k - 1) % m == 0 { out.push((-6 + ((k - 1) / m) as i32, y.clone())); } // t = 10^e
        }
        println!("backward Euler, {} steps per decade ({} steps, {:.3} s):", m, times.len() - 1,
                 t.elapsed().as_secs_f64());
        for (e, yk) in out.iter().filter(|(e, _)| [-3, 0, 4, 11].contains(e)) {
            println!("  1e{}  {:.4e}  {:.4e}  {:.4}   sum - 1 = {:.1e}", e, yk[0], yk[1], yk[2], yk.sum() - 1.0);
        }
    }
}
