// Q4: a double root has no sign change, so Brent's method (argmin) refuses to start.
use argmin::core::{CostFunction, Error, Executor};
use argmin::solver::brent::BrentRoot;

struct Square;
impl CostFunction for Square {
    type Param = f64;
    type Output = f64;
    fn cost(&self, x: &f64) -> Result<f64, Error> {
        Ok((x - 2.0).powi(2)) // (x - 2)^2: root at 2, but never negative
    }
}

fn main() {
    match Executor::new(Square, BrentRoot::new(0.0, 5.0, 1e-10)).run() {
        Ok(res) => println!("{}", res),
        Err(e) => println!("{e}"),
    }
}
