// C4(a): fit the logistic curve with Nelder-Mead from (r, K, N0) = (0.5, 1, 0.05).
// R's optim and scipy's minimize are replaced by the argmin crate's NelderMead.
// argmin takes different steps than R or SciPy, so it may stop at a different
// (also wrong) point.
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

fn main() -> Result<(), Error> {
    let t_min: Vec<f64> = (0..=12).map(|i| 120.0 * i as f64).collect(); // minutes
    let od = vec![0.001, 0.047, 0.139, 0.282, 0.581, 0.859, 1.049,
                  1.161, 1.152, 1.217, 1.182, 1.176, 1.185];
    let start = vec![0.5, 1.0, 0.05];
    // argmin needs the whole starting simplex; build it like R's optim does:
    // the start, plus the start with 0.1 added to one parameter at a time
    let mut simplex = vec![start.clone()];
    for i in 0..3 {
        let mut v = start.clone();
        v[i] += 0.1;
        simplex.push(v);
    }
    // stop when the SSE values at the simplex corners differ by less than 1e-10
    let solver = NelderMead::new(simplex).with_sd_tolerance(1e-10)?;
    let res = Executor::new(Growth { t: t_min, od }, solver)
        .configure(|state| state.max_iters(500))
        .run()?; // ? stops main and returns the error if the fit failed
    let s = res.state();
    let p = s.get_best_param().unwrap();
    println!("{}  SSE {:.4}  r {:.4}  K {:.4}  N0 {:.3e}",
             s.get_termination_reason().unwrap().text(), s.get_best_cost(), p[0], p[1], p[2]);
    println!("doubling time {:.2} min", 2f64.ln() / p[0]);
    Ok(())
}
