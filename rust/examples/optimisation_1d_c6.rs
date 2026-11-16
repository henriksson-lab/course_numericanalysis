// C6: the growth rate r of the logistic OD curve, with K and N0 fixed (BrentOpt).
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::BrentOpt;

const K: f64 = 1.2;
const N0: f64 = 0.05;
const Y: [f64; 61] = [0.047, 0.058, 0.060, 0.081, 0.084, 0.090, 0.108, 0.123, 0.137, 0.150, 0.177,
    0.192, 0.210, 0.226, 0.269, 0.292, 0.322, 0.359, 0.393, 0.428, 0.466, 0.504,
    0.540, 0.569, 0.622, 0.659, 0.698, 0.730, 0.772, 0.812, 0.851, 0.876, 0.909,
    0.935, 0.955, 0.985, 1.007, 1.029, 1.054, 1.069, 1.079, 1.092, 1.109, 1.120,
    1.123, 1.132, 1.145, 1.154, 1.155, 1.166, 1.168, 1.167, 1.175, 1.171, 1.187,
    1.192, 1.183, 1.181, 1.191, 1.189, 1.203];

fn od(t: f64, r: f64) -> f64 { K / (1.0 + (K / N0 - 1.0) * (-r * t).exp()) }
fn sse(r: f64) -> f64 {
    Y.iter().enumerate().map(|(i, y)| (y - od(i as f64 / 6.0, r)).powi(2)).sum() // t = i/6 h
}

struct Sse;
impl CostFunction for Sse {
    type Param = f64;
    type Output = f64;
    fn cost(&self, r: &f64) -> Result<f64, Error> { Ok(sse(*r)) }
}

fn main() -> Result<(), Error> {
    let look: Vec<String> = [0.1, 0.4, 0.8, 1.2, 2.0, 3.0].iter().map(|&r| format!("{:.3}", sse(r))).collect();
    println!("{}", look.join(" ")); // (a) look first
    let solver = BrentOpt::new(0.1, 3.0).set_tolerance(f64::EPSILON.sqrt(), 1e-10); // (b)
    let res = Executor::new(Sse, solver).configure(|s| s.max_iters(200)).run()?;
    let r = *res.state().get_best_param().unwrap();
    println!("{:.9} {:.9}", r, sse(r));
    println!("{:.9}", (sse(r) / 60.0).sqrt()); // (c) residual sd
    Ok(())
}
