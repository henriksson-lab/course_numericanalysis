// rust/examples/iterative_solvers_c1.rs  —  cargo run --release --example iterative_solvers_c1
// C1: Jacobi and Gauss-Seidel as splittings x <- x + M^{-1} (b - A x).
// Gauss-Seidel's M = lower triangle of A, solved by forward substitution.
include!("iterative_solvers_shared/defs.rs");

fn splitting(a: &CsMat<f64>, b: &[f64], msolve: impl Fn(&[f64]) -> Vec<f64>,
             tol: f64, maxit: usize) -> (Vec<f64>, usize, Vec<f64>) {
    let mut x = vec![0.0; b.len()];
    let (nb, mut res) = (norm(b), Vec::new());
    let mut it = 0;
    loop {
        let ax = matvec(a, &x);
        let r: Vec<f64> = b.iter().zip(&ax).map(|(p, q)| p - q).collect();
        res.push(norm(&r) / nb);
        if *res.last().unwrap() < tol || it == maxit { break; }
        let dx = msolve(&r);
        for i in 0..x.len() { x[i] += dx[i]; }
        it += 1;
    }
    (x, it, res)
}
fn forward_solve(a: &CsMat<f64>, r: &[f64]) -> Vec<f64> {   // (lower triangle of A) y = r
    let mut y = vec![0.0; r.len()];
    for (i, row) in a.outer_iterator().enumerate() {
        let (mut s, mut d) = (r[i], 0.0);
        for (j, v) in row.iter() {
            if j < i { s -= v * y[j]; } else if j == i { d = *v; }
        }
        y[i] = s / d;
    }
    y
}

fn main() {
    let mut rj = Vec::new();
    for m in [7, 15, 31] {
        let (a, b) = (lap2(m), chamber_rhs(m, 12.0));
        let h = 1.0 / (m + 1) as f64;
        let d = diagonal(&a);
        let (_, itj, r) = splitting(&a, &b, |r| r.iter().zip(&d).map(|(p, q)| p / q).collect(), 1e-6, 100_000);
        let (_, itg, _) = splitting(&a, &b, |r| forward_solve(&a, r), 1e-6, 100_000);
        let pred = 1e-6_f64.ln() / (std::f64::consts::PI * h).cos().ln();
        println!("{} {} {} {} {}", m, itj, itg, pred.round(), (pred / 2.0).round());
        rj = r;
    }
    let k = rj.len();
    let ratios: Vec<f64> = (k - 3..k).map(|i| (rj[i] / rj[i - 1] * 1e5).round() / 1e5).collect();
    println!("{:?}", ratios);                         // residual ratios, m = 31
    println!("{}", (std::f64::consts::PI / 32.0).cos());
}
