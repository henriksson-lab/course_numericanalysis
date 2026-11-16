// C5: the period of a circadian signal. RSS(T) from a linear least-squares fit
// (nalgebra SVD) for each period T; BrentOpt on a wide interval, then grid + refine.
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::BrentOpt;
use nalgebra::{DMatrix, DVector};
use std::f64::consts::PI;

const Y: [f64; 36] = [2.10, 3.07, 2.42, 2.29, 1.64, 1.08, 1.95, 2.81, 2.52, 1.77, 1.80, 1.32,
                      1.92, 2.06, 3.02, 2.17, 1.59, 1.42, 1.44, 2.22, 3.01, 2.67, 1.54, 1.27,
                      2.03, 2.44, 2.82, 2.11, 1.82, 1.67, 1.52, 2.24, 3.03, 2.60, 2.08, 1.55];

// best sine wave with period T: least squares on the columns 1, cos, sin
fn rss(period: f64) -> f64 {
    let n = Y.len();
    let x = DMatrix::from_fn(n, 3, |i, j| {
        let t = 4.0 * i as f64; // hours: 0, 4, ..., 140
        match j { 0 => 1.0, 1 => (2.0 * PI * t / period).cos(), _ => (2.0 * PI * t / period).sin() }
    });
    let y = DVector::from_column_slice(&Y);
    let coef = x.clone().svd(true, true).solve(&y, 1e-12).unwrap();
    (y - x * coef).norm_squared()
}

struct Rss;
impl CostFunction for Rss {
    type Param = f64;
    type Output = f64;
    fn cost(&self, p: &f64) -> Result<f64, Error> { Ok(rss(*p)) }
}

fn brent(a: f64, b: f64, t: f64) -> Result<f64, Error> {
    let solver = BrentOpt::new(a, b).set_tolerance(f64::EPSILON.sqrt(), t);
    let res = Executor::new(Rss, solver).configure(|s| s.max_iters(200)).run()?;
    Ok(*res.state().get_best_param().unwrap())
}

fn main() -> Result<(), Error> {
    println!("{:.5} {:.5}", brent(10.0, 40.0, 1e-5)?, brent(15.0, 35.0, 1e-5)?);

    // (b) grid, then refine
    let grid: Vec<f64> = (0..=120).map(|i| 10.0 + 0.25 * i as f64).collect();
    let t0 = grid.iter().cloned()
        .min_by(|a, b| rss(*a).partial_cmp(&rss(*b)).unwrap()).unwrap();
    println!("{}", t0);
    let tbest = brent(t0 - 0.5, t0 + 0.5, 1e-6)?;
    println!("{:.5} {:.6} {:.6}", tbest, rss(tbest), rss(32.40162));

    // (c) a coarse grid
    let coarse: Vec<String> = [10.0, 20.0, 30.0, 40.0].iter().map(|&p| format!("{:.2}", rss(p))).collect();
    println!("{}", coarse.join(" "));
    Ok(())
}
