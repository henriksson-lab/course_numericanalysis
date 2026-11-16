// C4(c): the same fit with time in hours instead of minutes.
// R's optim and scipy's minimize are replaced by the argmin crate's NelderMead.
// argmin takes different steps than R or SciPy, so the last digits (and, from
// a bad start, where it gets stuck) can differ.
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::neldermead::NelderMead;

fn logistic(p: &[f64], t: f64) -> f64 {
    let (r, k, n0) = (p[0], p[1], p[2]); // p = (r, K, N0)
    k / (1.0 + (k / n0 - 1.0) * (-r * t).exp())
}

// argmin wants the problem as a type with a `cost` method: here, the SSE
struct Growth { t: Vec<f64>, od: Vec<f64> }

impl CostFunction for Growth {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Vec<f64>) -> Result<f64, Error> {
        Ok(self.t.iter().zip(&self.od).map(|(&t, &y)| (y - logistic(p, t)).powi(2)).sum())
    }
}

// Nelder-Mead from `start`, with R's optim starting simplex: the start, plus
// the start with 0.1 added to one parameter at a time
fn fit(start: &[f64], t: &[f64], od: &[f64]) -> Result<(String, f64, Vec<f64>), Error> {
    let mut simplex = vec![start.to_vec()];
    for i in 0..start.len() {
        let mut v = start.to_vec();
        v[i] += 0.1;
        simplex.push(v);
    }
    // stop when the SSE values at the simplex corners differ by less than 1e-10
    let solver = NelderMead::new(simplex).with_sd_tolerance(1e-10)?;
    let problem = Growth { t: t.to_vec(), od: od.to_vec() };
    let res = Executor::new(problem, solver).configure(|s| s.max_iters(500)).run()?;
    let s = res.state();
    let reason = s.get_termination_reason().unwrap().text().to_string();
    Ok((reason, s.get_best_cost(), s.get_best_param().unwrap().clone()))
}

fn main() -> Result<(), Error> {
    let od = [0.001, 0.047, 0.139, 0.282, 0.581, 0.859, 1.049,
              1.161, 1.152, 1.217, 1.182, 1.176, 1.185];
    let t_h: Vec<f64> = (0..=12).map(|i| 120.0 * i as f64 / 60.0).collect(); // hours
    let (reason, sse, p) = fit(&[0.5, 1.0, 0.05], &t_h, &od)?;
    println!("{}  SSE {:.4}  r {:.4}  K {:.3}  N0 {:.5}", reason, sse, p[0], p[1], p[2]);
    println!("doubling time {:.1} min", 60.0 * 2f64.ln() / p[0]);
    Ok(())
}
