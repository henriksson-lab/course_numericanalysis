// C6: heat conduction through the wall of a PCR chip: PDMS vs glass, 1 vs 2 mm.
// Backward Euler with Thomas; heater at z = 0, channel side at z = H
// (insulated, or losing heat to air through a Robin condition).

/// Thomas algorithm (a[0] and c[n-1] unused), as in Structure: symmetric, banded, sparse.
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec());
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    x
}

/// Top (channel-side) temperature after every step. Units: mm, s, W, K.
fn wall(h_mm: f64, a: f64, k: f64, h: f64, tend: f64) -> Vec<f64> {
    let (n, dt, tair) = (50usize, 0.05, 25.0);
    let dz = h_mm / n as f64;
    let r = a * dt / (dz * dz);
    let m = n;                                           // unknowns z_1 .. z_N
    let beta = 2.0 * dz * h / (k * 1000.0);              // Robin ghost point (dz in mm)
    let mut lo = vec![-r; m]; lo[0] = 0.0; lo[m - 1] = -2.0 * r;
    let mut di = vec![1.0 + 2.0 * r; m]; di[m - 1] += r * beta;
    let mut up = vec![-r; m]; up[m - 1] = 0.0;
    let mut t = vec![65.0; m];
    let mut top = Vec::new();
    for _ in 0..(tend / dt).round() as usize {           // heater at 95 C from t = 0
        let mut d = t.clone();
        d[0] += r * 95.0;
        d[m - 1] += r * beta * tair;
        t = thomas(&lo, &di, &up, &d);
        top.push(t[m - 1]);
    }
    top
}

fn main() {
    let mats = [("PDMS", 0.15, 970.0, 1460.0), ("glass", 1.1, 2200.0, 700.0)];
    for (name, k, rho, cp) in mats {
        let a = k / (rho * cp) * 1e6;                    // mm^2/s
        let reach: Vec<f64> = [1.0, 2.0].iter().map(|&hh| {
            let top = wall(hh, a, k, 0.0, 120.0);
            0.05 * (top.iter().position(|&v| v >= 94.5).unwrap() + 1) as f64
        }).collect();
        let drop: Vec<f64> = [1.0, 2.0].iter()
            .map(|&hh| 95.0 - wall(hh, a, k, 10.0, 600.0).last().unwrap()).collect();
        println!("{:5}  alpha {:.4} mm^2/s  reach 94.5 C: {:.2} s, {:.2} s  steady drop: {:.3} C, {:.3} C",
                 name, a, reach[0], reach[1], drop[0], drop[1]);
    }
}
