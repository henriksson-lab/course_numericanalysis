// rust/examples/gradient_methods_c2.rs  —  cargo run --release --example gradient_methods_c2
// C2: steepest descent with Armijo backtracking on f = (x1^2 + kappa x2^2)/2.
fn sd_armijo(f: &dyn Fn(&[f64]) -> f64, g: &dyn Fn(&[f64]) -> Vec<f64>,
             mut x: Vec<f64>, tol: f64, maxit: usize) -> (Vec<f64>, usize) {
    let mut k = 0;
    while k < maxit {
        let gk = g(&x);
        let g2: f64 = gk.iter().map(|v| v * v).sum();
        if g2.sqrt() < tol { break; }
        let mut a = 1.0;
        let fx = f(&x);
        let trial = |a: f64| -> Vec<f64> { x.iter().zip(&gk).map(|(x, g)| x - a * g).collect() };
        // `!(.. <= ..)` also rejects NaN
        while !(f(&trial(a)) <= fx - 1e-4 * a * g2) { a /= 2.0; }
        x = trial(a);
        k += 1;
    }
    (x, k)
}

fn main() {
    for kap in [1.0, 10.0, 100.0, 1000.0f64] {
        let f = move |x: &[f64]| (x[0] * x[0] + kap * x[1] * x[1]) / 2.0;
        let g = move |x: &[f64]| vec![x[0], kap * x[1]];
        let (_, steps) = sd_armijo(&f, &g, vec![1.0, 1.0], 1e-8, 100_000);
        let pred = if kap == 1.0 { 1.0 } else { 1e-8f64.ln() / ((kap - 1.0) / (kap + 1.0)).ln() };
        println!("kappa {:5}  steps {:6}  predicted {:8.0}", kap, steps, pred);
    }
}
