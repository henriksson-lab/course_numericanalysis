// C6: estimate D by fitting the marching model to a synthetic T-sensor profile, with
// Brent's method (argmin's BrentOpt). Grid refinement of the fit; a wrong flow rate.
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::BrentOpt;

fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec());
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    x
}

// backward Euler in x (as in C3)
fn march(c: &[f64], u: &[f64], d: f64, dy: f64, dx: f64, nsteps: usize) -> Vec<f64> {
    let m = c.len();
    let a: Vec<f64> = u.iter().map(|uj| d * dx / (uj * dy * dy)).collect();
    let mut di: Vec<f64> = a.iter().map(|aj| 1.0 + 2.0 * aj).collect();
    di[0] = 1.0 + a[0];
    di[m - 1] = 1.0 + a[m - 1];
    let off: Vec<f64> = a.iter().map(|aj| -aj).collect();
    let mut c = c.to_vec();
    for _ in 0..nsteps {
        c = thomas(&off, &di, &off, &c);
    }
    c
}

const COBS: [f64; 25] = [0.881, 0.906, 0.864, 0.859, 0.852, 0.828, 0.824, 0.758, 0.753, 0.681, 0.617,
                         0.578, 0.496, 0.457, 0.383, 0.346, 0.269, 0.216, 0.153, 0.157, 0.112, 0.180,
                         0.100, 0.089, 0.125];
const W: f64 = 100.0;

// model at x = 1 mm, averaged over 25 pixels of 4 um
fn profile(d: f64, m: usize, dx: f64, ubar: f64) -> Vec<f64> {
    let dy = W / m as f64;
    let y: Vec<f64> = (1..=m).map(|j| (j as f64 - 0.5) * dy).collect();
    let u: Vec<f64> = y.iter().map(|yj| 6.0 * ubar * (yj / W) * (1.0 - yj / W)).collect();
    let c0: Vec<f64> = y.iter().map(|&yj| if yj < W / 2.0 { 1.0 } else { 0.0 }).collect();
    let c = march(&c0, &u, d, dy, dx, (1000.0 / dx).round() as usize);
    c.chunks(m / 25).map(|p| p.iter().sum::<f64>() / p.len() as f64).collect()
}

struct Sse { m: usize, dx: f64, ubar: f64 }
impl CostFunction for Sse {
    type Param = f64;
    type Output = f64;
    fn cost(&self, d: &f64) -> Result<f64, Error> {
        Ok(profile(*d, self.m, self.dx, self.ubar).iter().zip(COBS.iter()).map(|(p, o)| (o - p).powi(2)).sum())
    }
}

fn fit_d(m: usize, dx: f64, ubar: f64) -> Result<(f64, f64), Error> {
    let solver = BrentOpt::new(10.0, 3000.0).set_tolerance(f64::EPSILON.sqrt(), 1e-6);
    let res = Executor::new(Sse { m, dx, ubar }, solver).configure(|s| s.max_iters(100)).run()?;
    Ok((*res.state().get_best_param().unwrap(), res.state().get_best_cost()))
}

fn main() -> Result<(), Error> {
    let (d, s) = fit_d(50, 25.0, 1000.0)?;
    println!("fit: D = {:.2} um^2/s, SSE = {:.5}", d, s);
    for (m, dx) in [(50, 1000.0), (50, 250.0), (50, 100.0), (50, 25.0), (100, 10.0), (200, 5.0)] {
        println!("M = {:3}, dx = {:6.1}: D = {:.2}", m, dx, fit_d(m, dx, 1000.0)?.0);
    }
    println!("if the flow were 800 um/s: D = {:.2}", fit_d(50, 25.0, 800.0)?.0);
    Ok(())
}
