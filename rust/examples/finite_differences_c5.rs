// C5: local least-squares slope over a window of 2m+1 points.
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

// slope at points m..n-m-1: sum_j j*y[i+j] / (dt * sum_j j^2)
fn local_slope(y: &[f64], m: usize, dt: f64) -> Vec<f64> {
    let m = m as i64;
    let sjj: f64 = (-m..=m).map(|j| (j * j) as f64).sum();
    (m..y.len() as i64 - m)
        .map(|i| (-m..=m).map(|j| j as f64 * y[(i + j) as usize]).sum::<f64>() / (dt * sjj))
        .collect()
}

fn main() {
    let dt = 1.0 / 6.0;
    let mut rms = Vec::new();
    for m in 1..=8 {
        let d = local_slope(&Y, m, dt);
        let ss: f64 = d.iter().enumerate()
            .map(|(k, d)| (d - gr((k + m) as f64 * dt)).powi(2)).sum();
        rms.push((ss / d.len() as f64).sqrt());
        println!("{} {:.3e}", m, rms[m - 1]);
    }
    let m = 1 + (0..8).min_by(|&a, &b| rms[a].partial_cmp(&rms[b]).unwrap()).unwrap();
    let d = local_slope(&Y, m, dt);
    let imax = (0..d.len()).max_by(|&a, &b| d[a].partial_cmp(&d[b]).unwrap()).unwrap();
    println!("m {}  max {:.4} at {:.4} h", m, d[imax], (imax + m) as f64 * dt);
}
