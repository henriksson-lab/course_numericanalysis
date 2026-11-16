// rust/examples/iterative_solvers_c3.rs  —  cargo run --release --example iterative_solvers_c3
// C3: hand-written CG (pcg with no preconditioner, in the shared defs) on sprs:
// iterations vs grid size, and the slope on log-log axes.
include!("iterative_solvers_shared/defs.rs");
use std::f64::consts::PI;

fn main() {
    let ms = [31usize, 63, 127, 255];
    let mut its = Vec::new();
    for &m in &ms {
        let (a, b) = (lap2(m), chamber_rhs(m, 12.0));
        let (x, it) = pcg(&a, &b, |r| r.to_vec(), 1e-6, 10_000);
        let it = it.unwrap();
        its.push(it as f64);
        let kappa = 1.0 / (PI / (2.0 * (m + 1) as f64)).tan().powi(2);   // cot^2(pi h / 2)
        println!("{} {} {} {} {}", m, it, kappa.round(), kappa.sqrt().round(), 1.0 - x[(m * m - 1) / 2]);
    }
    // least-squares slope of log(its) against log(m + 1)
    let lx: Vec<f64> = ms.iter().map(|&m| ((m + 1) as f64).ln()).collect();
    let ly: Vec<f64> = its.iter().map(|v| v.ln()).collect();
    let (mx, my) = (lx.iter().sum::<f64>() / 4.0, ly.iter().sum::<f64>() / 4.0);
    let sxy: f64 = lx.iter().zip(&ly).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = lx.iter().map(|x| (x - mx).powi(2)).sum();
    println!("{}", sxy / sxx);
}
