// Background: central difference at the best step, and a Richardson table.
// Rust has no numDeriv / scipy.differentiate in our crate list: both are hand-written.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}
fn gr(t: f64) -> f64 {
    R * od(t) * (1.0 - od(t) / K) // exact N'(t)
}

fn main() {
    let h = f64::EPSILON.powf(1.0 / 3.0); // best step for central
    let cen = (od(3.0 + h) - od(3.0 - h)) / (2.0 * h);
    println!("{} {:e}", cen, cen - gr(3.0));
    // Richardson, 4 levels from h = 1e-4 * |x| (what numDeriv::grad does)
    let mut d: Vec<f64> = (0..4)
        .map(|i| { let h = 3e-4 / 2f64.powi(i); (od(3.0 + h) - od(3.0 - h)) / (2.0 * h) })
        .collect();
    for j in 1..4 {
        let f = 4f64.powi(j as i32);
        for i in (j..4).rev() { d[i] = (f * d[i] - d[i - 1]) / (f - 1.0); }
    }
    println!("{:e}", d[3] - gr(3.0));
}
