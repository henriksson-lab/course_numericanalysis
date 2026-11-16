// C7: Brent's method from the argmin crate (Rust's counterpart of uniroot / brentq).
use argmin::core::{CostFunction, Error, Executor, State};
use argmin::solver::brent::BrentRoot;

const R: f64 = 0.5;
const K: f64 = 100.0;

// argmin wants the function as a struct implementing CostFunction
struct Steady { h: f64 }
impl CostFunction for Steady {
    type Param = f64;
    type Output = f64;
    fn cost(&self, n: &f64) -> Result<f64, Error> {
        Ok(n * (R * (1.0 - n / K)).exp() - self.h - n) // f(N) = g(N) - N
    }
}

fn brent(h: f64, a: f64, b: f64, tol: f64) -> Result<(f64, u64), Error> {
    let res = Executor::new(Steady { h }, BrentRoot::new(a, b, tol))
        .configure(|state| state.max_iters(100))
        .run()?;                               // `?`: pass an error on to the caller
    let st = res.state();
    Ok((*st.get_best_param().unwrap(), st.get_iter()))
}

fn main() {
    for (a, b) in [(0.0, 50.0), (50.0, 100.0)] {
        let (x1, it1) = brent(10.0, a, b, 1e-4).unwrap();
        let (x2, it2) = brent(10.0, a, b, 1e-14).unwrap();
        println!("[{a}, {b}]: tol 1e-4: {x1:.10} ({it1} iter); tol 1e-14: {x2:.10} ({it2} iter)");
    }
    match brent(15.0, 0.0, 100.0, 1e-14) {
        Ok(_) => println!("found a root"),
        Err(e) => println!("{e}"),
    }
}
