// Background: degree-4 interpolant through five OD readings, and piecewise linear.
// Rust has no interpolation in std or our crates: both are hand-written.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

// barycentric formula (see the last lecture slide)
fn interp(x: &[f64], y: &[f64], z: f64) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for j in 0..x.len() {
        if z == x[j] {
            return y[j];
        }
        let mut w = 1.0;
        for k in 0..x.len() {
            if k != j {
                w /= x[j] - x[k];
            }
        }
        num += w * y[j] / (z - x[j]);
        den += w / (z - x[j]);
    }
    num / den
}

// piecewise linear, like R's approx() / np.interp (x increasing)
fn linear(x: &[f64], y: &[f64], z: f64) -> f64 {
    let i = (1..x.len()).find(|&i| z <= x[i]).unwrap();
    y[i - 1] + (y[i] - y[i - 1]) * (z - x[i - 1]) / (x[i] - x[i - 1])
}

fn main() {
    let t = [0.0, 2.0, 4.0, 6.0, 8.0]; // readings every 2 h
    let y: Vec<f64> = t.iter().map(|&ti| od(ti)).collect();
    for z in [3.0, 5.0] {
        println!("t = {z}: interpolant {:.7}  linear {:.7}  truth {:.7}",
                 interp(&t, &y, z), linear(&t, &y, z), od(z));
    }
}
