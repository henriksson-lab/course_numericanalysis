// C3: is glucose limiting at the closed end of a mother-machine trench in rich
// medium (0.2 % glucose = 11.1 mM)? Vary the effective diffusion coefficient.
// thomas and trench_mm are the same as in steady_state_2d_c2.rs.

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
    for d in [600.0_f64, 60.0, 6.0] {
        let (c, _) = trench_mm(11.1, 25.0, d, 0.6, 0.005, 100);
        println!("D = {:3} um^2/s: closed end {:.3} mM, Lmax = {:.0} um",
                 d, c[99], (2.0 * d * 11.1 / 0.6).sqrt());
    }
}
