// rust/examples/phase_plane_c6.rs  —  cargo run --example phase_plane_c6
// C6: near the bifurcation at alpha = 2: the slowest eigenvalue, how long a simulation
// needs, and what a finite-difference Jacobian says when the true eigenvalue is 0.
use nalgebra::Matrix2;
use ode_solvers::{Dopri5, System, Vector2};

struct Toggle { alpha: f64 }
impl System<f64, Vector2<f64>> for Toggle {
    fn system(&self, _t: f64, x: &Vector2<f64>, dx: &mut Vector2<f64>) {
        dx[0] = self.alpha / (1.0 + x[1] * x[1]) - x[0];
        dx[1] = self.alpha / (1.0 + x[0] * x[0]) - x[1];
    }
}

fn f(alpha: f64, x: [f64; 2]) -> [f64; 2] {
    [alpha / (1.0 + x[1] * x[1]) - x[0], alpha / (1.0 + x[0] * x[0]) - x[1]]
}

fn jac(alpha: f64, u: f64, v: f64) -> Matrix2<f64> {
    Matrix2::new(-1.0, -2.0 * alpha * v / (1.0 + v * v).powi(2),
                 -2.0 * alpha * u / (1.0 + u * u).powi(2), -1.0)
}

fn fd_jacobian(alpha: f64, x: [f64; 2]) -> Matrix2<f64> {
    let h = f64::EPSILON.sqrt();
    let fx = f(alpha, x);
    let mut jm = Matrix2::zeros();
    for j in 0..2 {
        let hj = h * x[j].abs().max(1.0);
        let mut xh = x;
        xh[j] += hj;
        let fh = f(alpha, xh);
        for i in 0..2 { jm[(i, j)] = (fh[i] - fx[i]) / hj; }
    }
    jm
}

fn lmax(m: &Matrix2<f64>) -> f64 {
    m.complex_eigenvalues().iter().map(|z| z.re).fold(f64::NEG_INFINITY, f64::max)
}

fn main() {
    for alpha in [3.0, 2.5, 2.1, 2.01, 2.001] {
        let r = (alpha * alpha - 4.0_f64).sqrt();
        let (uon, von) = ((alpha + r) / 2.0, (alpha - r) / 2.0); // "U on", exact
        let l = lmax(&jac(alpha, uon, von));
        let mut s = Dopri5::new(Toggle { alpha }, 0.0, 100.0, 1.0, Vector2::new(3.0, 0.5), 1e-10, 1e-12);
        s.integrate().expect("integration failed");
        let x = s.y_out().last().unwrap();
        let err = ((x[0] - uon).powi(2) + (x[1] - von).powi(2)).sqrt();
        let e = if err < 1e-8 { "< 1e-08".to_string() } else { format!("{:.1e}", err) };
        println!("alpha = {:.3}  lambda_max = {:8.5}  t_1% = {:7.1}  error at t = 100: {}",
                 alpha, l, 100f64.ln() / l.abs(), e);
    }
    for x in [[1.0, 1.0], [1.0 + 1e-6, 1.0 - 1e-6], [1.0 - 1e-6, 1.0 + 1e-6]] {
        println!("at ({:.6}, {:.6}): lambda_max from FD Jacobian = {:.2e}", x[0], x[1],
                 lmax(&fd_jacobian(2.0, x)));
    }
}
