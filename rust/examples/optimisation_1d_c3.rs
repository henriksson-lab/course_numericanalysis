// C3: the flat bottom. Root of the score (BrentRoot) versus minimum of minus the
// log-likelihood (BrentOpt), for the zero-truncated Poisson UMI counts.
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::{BrentOpt, BrentRoot};

const N: f64 = 100.0; // cells
const SX: f64 = 173.0; // total UMIs

fn nll(l: f64) -> f64 { -(SX * l.ln() - N * l - N * (-(-l).exp_m1()).ln()) } // minus log-lik
fn score(l: f64) -> f64 { SX / l - N / (-(-l).exp_m1()) }                    // = -nll'
fn dscore(l: f64) -> f64 { -SX / l.powi(2) + N * (-l).exp() / (-l).exp_m1().powi(2) }

struct Nll;
impl CostFunction for Nll {
    type Param = f64;
    type Output = f64;
    fn cost(&self, l: &f64) -> Result<f64, Error> { Ok(nll(*l)) }
}
struct Score;
impl CostFunction for Score {
    type Param = f64;
    type Output = f64;
    fn cost(&self, l: &f64) -> Result<f64, Error> { Ok(score(*l)) }
}

fn main() -> Result<(), Error> {
    let r = Executor::new(Score, BrentRoot::new(0.01, 10.0, 1e-15))
        .configure(|s| s.max_iters(200)).run()?;
    let lhat = *r.state().get_best_param().unwrap();
    let o = Executor::new(Nll, BrentOpt::new(0.01, 10.0).set_tolerance(1e-15, 1e-15))
        .configure(|s| s.max_iters(200)).run()?;
    let lopt = *o.state().get_best_param().unwrap();
    println!("{:.16} {:.16}", lhat, lopt);
    println!("{:.6e}", lopt - lhat);
    for k in 4..=10 {
        let d = 10f64.powi(-k);
        println!("{:.0e}  {:>13.6e}  {:>13.6e}", d, nll(lhat + d) - nll(lhat), score(lhat + d));
    }
    // predicted size of the flat bottom: sqrt(2 eps |f| / f'')
    println!("{:.6e}", (2.0 * f64::EPSILON * nll(lhat) / -dscore(lhat)).sqrt());
    Ok(())
}
