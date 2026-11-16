// C3: the ATP LP with an oxygen limit (minilp). minilp has no shadow prices,
// so each one is measured by re-solving with a slightly larger right-hand side.
use minilp::{ComparisonOp, OptimizationDirection, Problem};

fn solve(rhs: [f64; 3]) -> (f64, f64, f64) {
    let mut p = Problem::new(OptimizationDirection::Maximize);
    let x1 = p.add_var(6.0, (0.0, f64::INFINITY));
    let x2 = p.add_var(3.0, (0.0, f64::INFINITY));
    p.add_constraint(&[(x1, 1.0), (x2, 1.0)], ComparisonOp::Le, rhs[0]); // carbon
    p.add_constraint(&[(x1, 4.0), (x2, 1.0)], ComparisonOp::Le, rhs[1]); // proteome
    p.add_constraint(&[(x1, 1.0)], ComparisonOp::Le, rhs[2]); // oxygen
    let sol = p.solve().unwrap();
    (sol[x1], sol[x2], sol.objective())
}

fn main() {
    let rhs = [6.0, 12.0, 1.5];
    let (x1, x2, z) = solve(rhs);
    println!("{} {}  ATP {}", x1, x2, z);
    let h = 0.01;
    let prices: Vec<String> = (0..3).map(|i| {
        let mut r = rhs;
        r[i] += h;
        format!("{:.3}", (solve(r).2 - z) / h)
    }).collect();
    println!("shadow prices (carbon, proteome, oxygen): {}", prices.join(" "));
}
