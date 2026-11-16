// C4: growth rate from noisy plate-reader data, central differences of width k.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;
const Y: [f64; 61] = [
    0.047, 0.058, 0.060, 0.081, 0.084, 0.090, 0.108, 0.123, 0.137, 0.150,
    0.177, 0.192, 0.210, 0.226, 0.269, 0.292, 0.322, 0.359, 0.393, 0.428,
    0.466, 0.504, 0.540, 0.569, 0.622, 0.659, 0.698, 0.730, 0.772, 0.812,
    0.851, 0.876, 0.909, 0.935, 0.955, 0.985, 1.007, 1.029, 1.054, 1.069,
    1.079, 1.092, 1.109, 1.120, 1.123, 1.132, 1.145, 1.154, 1.155, 1.166,
    1.168, 1.167, 1.175, 1.171, 1.187, 1.192, 1.183, 1.181, 1.191, 1.189,
    1.203];

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}
fn gr(t: f64) -> f64 {
    R * od(t) * (1.0 - od(t) / K)
}

// central difference with half-width k*dt: (times, estimates)
fn wide(k: usize, dt: f64) -> (Vec<f64>, Vec<f64>) {
    (k..Y.len() - k)
        .map(|i| (i as f64 * dt, (Y[i + k] - Y[i - k]) / (2.0 * k as f64 * dt)))
        .unzip()
}
fn rms(t: &[f64], d: &[f64]) -> f64 {
    (t.iter().zip(d).map(|(t, d)| (d - gr(*t)).powi(2)).sum::<f64>() / t.len() as f64).sqrt()
}

fn main() {
    let dt = 1.0 / 6.0;
    let (t1, d1) = wide(1, dt);
    let min = d1.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = d1.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("rms {:.6}  min {:.3}  max {:.3}", rms(&t1, &d1), min, max);
    let r: Vec<f64> = (1..=8).map(|k| { let (t, d) = wide(k, dt); rms(&t, &d) }).collect();
    for k in 1..=8 { println!("{} {} {:.3e}", k, 10 * k, r[k - 1]); }
    let kbest = 1 + (0..8).min_by(|&a, &b| r[a].partial_cmp(&r[b]).unwrap()).unwrap();
    let (t, d) = wide(kbest, dt);
    let imax = (0..d.len()).max_by(|&a, &b| d[a].partial_cmp(&d[b]).unwrap()).unwrap();
    println!("max {:.5} at {:.4} h", d[imax], t[imax]);
}
