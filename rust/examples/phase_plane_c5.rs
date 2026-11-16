// rust/examples/phase_plane_c5.rs  —  cargo run --example phase_plane_c5
// C5: the Sel'kov glycolysis model: eigenvalues at the fixed point for several b,
// then amplitude and period of the limit cycle from simulations.
use nalgebra::Matrix2;
use ode_solvers::{Dopri5, System, Vector2};

const A: f64 = 0.08;

struct Selkov { b: f64 }
impl System<f64, Vector2<f64>> for Selkov {
    fn system(&self, _t: f64, x: &Vector2<f64>, dx: &mut Vector2<f64>) {
        let (u, v) = (x[0], x[1]);
        dx[0] = -u + A * v + u * u * v;
        dx[1] = self.b - A * v - u * u * v;
    }
}

fn jsel(u: f64, v: f64) -> Matrix2<f64> {
    Matrix2::new(-1.0 + 2.0 * u * v, A + u * u, -2.0 * u * v, -(A + u * u))
}

// amplitude and period after the transient (t > 100)
fn cycle(u0: f64, v0: f64, b: f64, tol: f64) {
    let mut s = Dopri5::new(Selkov { b }, 0.0, 200.0, 0.01, Vector2::new(u0, v0), tol, tol);
    s.integrate().expect("integration failed");
    let (ts, xs) = (s.x_out(), s.y_out());
    let (tt, u): (Vec<f64>, Vec<f64>) = ts.iter().zip(xs).filter(|(t, _)| **t > 100.0)
        .map(|(t, x)| (*t, x[0])).unzip();
    let mut tc = Vec::new(); // upward crossings of u = b, linear interpolation
    for i in 0..u.len() - 1 {
        if u[i] < b && u[i + 1] >= b {
            tc.push(tt[i] + (b - u[i]) / (u[i + 1] - u[i]) * 0.01);
        }
    }
    let (lo, hi) = u.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &x| (l.min(x), h.max(x)));
    let per = if tc.len() > 1 {
        format!("{:.4}", (tc[tc.len() - 1] - tc[0]) / (tc.len() - 1) as f64) // mean period
    } else { "none".to_string() };
    println!("({:.2}, {:.2}) b = {:.1} tol = {:.0e}: u in [{:.4}, {:.4}], range {:.1e}, period {}",
             u0, v0, b, tol, lo, hi, hi - lo, per);
}

fn main() {
    for b in [0.2, 0.4, 0.6, 0.8, 1.0] {
        let lam = jsel(b, b / (A + b * b)).complex_eigenvalues();
        println!("b = {:.1}  lambda = {:.4} +{:.4}i", b, lam[0].re, lam[0].im.abs());
    }
    cycle(0.61, 1.36, 0.6, 1e-6);
    cycle(2.5, 0.2, 0.6, 1e-6);
    cycle(2.5, 0.2, 0.6, 1e-10);
    cycle(2.5, 0.2, 1.0, 1e-10);
}
