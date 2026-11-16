// Background: the ATP linear program with the minilp crate (Rust's lpSolve::lp / linprog).
use minilp::{ComparisonOp, OptimizationDirection, Problem};

fn main() {
    let mut p = Problem::new(OptimizationDirection::Maximize);
    // add_var(objective coefficient, (lower bound, upper bound))
    let x1 = p.add_var(6.0, (0.0, f64::INFINITY)); // respiration flux
    let x2 = p.add_var(3.0, (0.0, f64::INFINITY)); // fermentation flux
    p.add_constraint(&[(x1, 1.0), (x2, 1.0)], ComparisonOp::Le, 6.0); // carbon
    p.add_constraint(&[(x1, 4.0), (x2, 1.0)], ComparisonOp::Le, 12.0); // proteome
    let sol = p.solve().unwrap(); // Err(Infeasible) or Err(Unbounded) if it fails
    println!("{} {}", sol[x1], sol[x2]);
    println!("{}", sol.objective());
    // minilp reports no shadow prices; C3 gets them by re-solving
}
