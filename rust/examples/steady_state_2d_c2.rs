// C2: nutrient (glucose) in a dead-end mother-machine trench with
// Michaelis-Menten uptake, solved by Newton. The Jacobian is tridiagonal,
// so each Newton step is one Thomas solve.

fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = b.len();
    let (mut bb, mut dd) = (b.to_vec(), d.to_vec());
    for i in 1..n {
        let w = a[i] / bb[i - 1];
        bb[i] -= w * c[i - 1];
        dd[i] -= w * dd[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = dd[n - 1] / bb[n - 1];
    for i in (0..n - 1).rev() { x[i] = (dd[i] - c[i] * x[i + 1]) / bb[i]; }
    x
}

// returns (concentrations at x_1..x_N, Newton steps)
fn trench_mm(c0: f64, l: f64, d: f64, vmax: f64, km: f64, n: usize) -> (Vec<f64>, usize) {
    let h = l / n as f64;
    let s = d / (h * h);
    let mut a = vec![-s; n];
    a[n - 1] = -2.0 * s;                                 // ghost point at the closed end
    let sup = vec![-s; n];
    let mut c = vec![c0; n];                             // start: trench = channel
    let mut it = 0;
    loop {
        // residual F = A c + r(c) - b; the known c0 at x_0 is the b part
        let f: Vec<f64> = (0..n).map(|i| {
            let left = if i == 0 { c0 } else { c[i - 1] };
            let right = if i + 1 < n { c[i + 1] } else { c[i - 1] };   // ghost point
            s * (-left + 2.0 * c[i] - right) + vmax * c[i] / (km + c[i])
        }).collect();
        let norm = f.iter().map(|v| v * v).sum::<f64>().sqrt();
        if norm < 1e-10 * (1.0 + s * c0) || it == 50 { break; }
        let diag: Vec<f64> = c.iter().map(|ci| 2.0 * s + vmax * km / ((km + ci) * (km + ci))).collect();
        let minus_f: Vec<f64> = f.iter().map(|v| -v).collect();
        let delta = thomas(&a, &diag, &sup, &minus_f);
        for i in 0..n { c[i] = (c[i] + delta[i]).max(0.0); }  // clip: no negative glucose
        it += 1;
    }
    (c, it)
}

fn main() {
    for c0 in [11.1, 1.0, 0.3, 0.1, 0.02] {
        let (c, it) = trench_mm(c0, 25.0, 600.0, 0.6, 0.005, 100);
        let starved = c.iter().filter(|&&v| v < 0.005).count() as f64 / 100.0 * 25.0;
        println!("c0 = {:5.2} mM: {} steps, closed end {:.4e} mM, starved {:4.1} um",
                 c0, it, c[99], starved);
    }
}
