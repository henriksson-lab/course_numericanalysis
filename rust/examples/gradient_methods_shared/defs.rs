// Shared by the gradient_methods_* examples: include!("gradient_methods_shared/defs.rs")
// The 4PL dose-response data, model, SSE and analytic gradient.
use argmin::core::{CostFunction, Error, Gradient};

#[allow(dead_code)]
const Y: [f64; 18] = [12.9, 9.8, 2.9, 0.1, 4.2, 9.7, 17.7, 22.3, 49.7, 42.1,
                      67.3, 72.6, 97.4, 91.8, 94.6, 93.3, 92.9, 95.1];
#[allow(dead_code)]
const P0: [f64; 4] = [0.0, 100.0, 1.0, 1.0];

#[allow(dead_code)]
fn doses() -> Vec<f64> {          // 9 doses 0.01 ... 100 µM, 2 wells each
    (0..18).map(|i| 10f64.powf(-2.0 + 0.5 * (i / 2) as f64)).collect()
}
#[allow(dead_code)]
fn f4(p: &[f64], c: f64) -> f64 {  // p = (bottom, top, EC50, hill)
    p[0] + (p[1] - p[0]) / (1.0 + (p[2] / c).powf(p[3]))
}
#[allow(dead_code)]
fn sse(p: &[f64]) -> f64 {
    doses().iter().zip(Y.iter()).map(|(&c, &y)| (y - f4(p, c)).powi(2)).sum()
}
#[allow(dead_code)]
fn grad_sse(p: &[f64]) -> Vec<f64> {   // -2 J^T r, one dose at a time
    let mut g = vec![0.0; 4];
    for (&c, &y) in doses().iter().zip(Y.iter()) {
        let u = (p[2] / c).powf(p[3]);
        let d = 1.0 + u;
        let r = y - f4(p, c);
        let j = [u / d, 1.0 / d,
                 -(p[1] - p[0]) * p[3] * u / (p[2] * d * d),
                 -(p[1] - p[0]) * u * (p[2] / c).ln() / (d * d)];
        for k in 0..4 { g[k] -= 2.0 * j[k] * r; }
    }
    g
}
#[allow(dead_code)]
fn norm(v: &[f64]) -> f64 { v.iter().map(|x| x * x).sum::<f64>().sqrt() }
#[allow(dead_code)]
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(x, y)| x * y).sum() }
// x + a * s
#[allow(dead_code)]
fn axpy(x: &[f64], a: f64, s: &[f64]) -> Vec<f64> { x.iter().zip(s).map(|(x, s)| x + a * s).collect() }

// argmin wants the problem as a type with a cost (and gradient) method
#[allow(dead_code)]
struct Fit;
impl CostFunction for Fit {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Vec<f64>) -> Result<f64, Error> { Ok(sse(p)) }
}
impl Gradient for Fit {
    type Param = Vec<f64>;
    type Gradient = Vec<f64>;
    fn gradient(&self, p: &Vec<f64>) -> Result<Vec<f64>, Error> { Ok(grad_sse(p)) }
}
