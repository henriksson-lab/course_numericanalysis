// Quiz Q9: a 2-variable LP with minilp; shadow prices by re-solving.
use minilp::{ComparisonOp, OptimizationDirection, Problem};

fn solve(b: [f64; 2]) -> (f64, f64, f64) {
    let mut p = Problem::new(OptimizationDirection::Maximize);
    let x1 = p.add_var(3.0, (0.0, f64::INFINITY));
    let x2 = p.add_var(2.0, (0.0, f64::INFINITY));
    p.add_constraint(&[(x1, 1.0), (x2, 1.0)], ComparisonOp::Le, b[0]);
    p.add_constraint(&[(x1, 1.0), (x2, 3.0)], ComparisonOp::Le, b[1]);
    let s = p.solve().unwrap();
    (s[x1], s[x2], s.objective())
}

fn main() {
    let (x1, x2, z) = solve([4.0, 6.0]);
    println!("{} {} {}", x1, x2, z);
    let h = 0.01;
    println!("{:.3} {:.3}", (solve([4.0 + h, 6.0]).2 - z) / h, (solve([4.0, 6.0 + h]).2 - z) / h);
}
