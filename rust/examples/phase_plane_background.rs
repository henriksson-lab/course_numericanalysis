// rust/examples/phase_plane_background.rs  —  cargo run --example phase_plane_background
// Background: one trajectory of the toggle switch (ode_solvers, Dopri5),
// and the eigenvalues of the Jacobian at "U on" (nalgebra).
use nalgebra::Matrix2;
use ode_solvers::{Dopri5, System, Vector2};

const ALPHA: f64 = 3.0;
const N: f64 = 2.0;

struct Toggle;
impl System<f64, Vector2<f64>> for Toggle {
    // ode_solvers wants system(t, x, dx): write the derivatives into dx
    fn system(&self, _t: f64, x: &Vector2<f64>, dx: &mut Vector2<f64>) {
        let (u, v) = (x[0], x[1]);
        dx[0] = ALPHA / (1.0 + v.powf(N)) - u; // du/dt
        dx[1] = ALPHA / (1.0 + u.powf(N)) - v; // dv/dt
    }
}

// Jacobian, as in Newton for systems
fn jac(u: f64, v: f64) -> Matrix2<f64> {
    Matrix2::new(-1.0, -ALPHA * N * v.powf(N - 1.0) / (1.0 + v.powf(N)).powi(2),
                 -ALPHA * N * u.powf(N - 1.0) / (1.0 + u.powf(N)).powi(2), -1.0)
}

fn main() {
    // Dopri5::new(system, t0, t_end, output spacing, x0, rtol, atol)
    let mut solver = Dopri5::new(Toggle, 0.0, 20.0, 1.0, Vector2::new(0.5, 0.3), 1e-10, 1e-12);
    solver.integrate().expect("integration failed");
    let (ts, xs) = (solver.x_out(), solver.y_out()); // output at t = 0, 1, 2, ..., 20
    for i in [0, 1, 5, 20] {
        println!("t = {:4.1}   u = {:.6}   v = {:.6}", ts[i], xs[i][0], xs[i][1]);
    }
    let (u_on, v_on) = ((3.0 + 5f64.sqrt()) / 2.0, (3.0 - 5f64.sqrt()) / 2.0); // "U on", exact
    println!("{:?}", jac(u_on, v_on).complex_eigenvalues().as_slice());
}
