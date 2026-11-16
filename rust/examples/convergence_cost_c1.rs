// C1: the V-shaped error of the forward difference.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

fn main() {
    let t0 = 3.0;
    let exact = R * od(t0) * (1.0 - od(t0) / K);

    let h: Vec<f64> = (1..=15).map(|k| 10f64.powi(-k)).collect();
    let err: Vec<f64> = h.iter().map(|&h| ((od(t0 + h) - od(t0)) / h - exact).abs()).collect();
    for (hi, e) in h.iter().zip(&err) {
        println!("{hi:.0e}  {e:.2e}"); // Rust writes 1e-1, not 1e-01
    }
    // index of the smallest error (total_cmp: a total order on f64)
    let best = (0..err.len()).min_by(|&i, &j| err[i].total_cmp(&err[j])).unwrap();
    println!("{:e}", h[best]);
    println!("{:e}", f64::EPSILON.sqrt());
    // No plotting here: write h and err to a file and plot in R or Python.
}
