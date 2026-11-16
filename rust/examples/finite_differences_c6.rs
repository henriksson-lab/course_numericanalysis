// C6: gradient of a log-likelihood by finite differences, against the exact one.
const SIGMA: f64 = 0.005;
const Y: [f64; 61] = [
    0.047, 0.058, 0.060, 0.081, 0.084, 0.090, 0.108, 0.123, 0.137, 0.150,
    0.177, 0.192, 0.210, 0.226, 0.269, 0.292, 0.322, 0.359, 0.393, 0.428,
    0.466, 0.504, 0.540, 0.569, 0.622, 0.659, 0.698, 0.730, 0.772, 0.812,
    0.851, 0.876, 0.909, 0.935, 0.955, 0.985, 1.007, 1.029, 1.054, 1.069,
    1.079, 1.092, 1.109, 1.120, 1.123, 1.132, 1.145, 1.154, 1.155, 1.166,
    1.168, 1.167, 1.175, 1.171, 1.187, 1.192, 1.183, 1.181, 1.191, 1.189,
    1.203];

fn t(i: usize) -> f64 { i as f64 / 6.0 }

fn loglik(th: &[f64]) -> f64 {
    let (r, k, n0) = (th[0], th[1], th[2]);
    -(0..61).map(|i| (Y[i] - k / (1.0 + (k / n0 - 1.0) * (-r * t(i)).exp())).powi(2))
        .sum::<f64>() / (2.0 * SIGMA * SIGMA)
}

fn grad_exact(th: &[f64], bug: bool) -> [f64; 3] {
    let (r, k, n0) = (th[0], th[1], th[2]);
    let mut g = [0.0; 3];
    for i in 0..61 {
        let e = (-r * t(i)).exp();
        let d = 1.0 + (k / n0 - 1.0) * e;
        let n = k / d;
        let w = (Y[i] - n) / (SIGMA * SIGMA);
        let tr = if bug { 1.0 } else { t(i) };
        g[0] += w * tr * n * (1.0 - n / k);            // dN/dr
        g[1] += w * (1.0 / d - k * e / (n0 * d * d));  // dN/dK
        g[2] += w * k * k * e / (n0 * n0 * d * d);     // dN/dN0
    }
    g
}

fn fd_grad(f: fn(&[f64]) -> f64, x: &[f64], h: &[f64], central: bool) -> [f64; 3] {
    let mut g = [0.0; 3];
    for j in 0..3 {
        let (mut xp, mut xm) = (x.to_vec(), x.to_vec());
        xp[j] += h[j];
        xm[j] -= h[j];
        g[j] = if central { (f(&xp) - f(&xm)) / (2.0 * h[j]) } else { (f(&xp) - f(x)) / h[j] };
    }
    g
}

fn main() {
    let th: [f64; 3] = [0.75, 1.15, 0.06];
    let g = grad_exact(&th, false);
    let rel = |a: [f64; 3]| -> Vec<String> {
        (0..3).map(|j| format!("{:.3e}", (a[j] - g[j]).abs() / g[j].abs())).collect()
    };
    let hrel = |c: f64| -> Vec<f64> { th.iter().map(|x| c * x.abs().max(1.0)).collect() };
    println!("fwd {:?}", rel(fd_grad(loglik, &th, &hrel(1e-8), false)));
    println!("cen {:?}", rel(fd_grad(loglik, &th, &hrel(6e-6), true)));
    // (b) at the maximum
    let mle: [f64; 3] = [0.79542677, 1.20245707, 0.050760766];
    let scaled = |c: f64| -> Vec<f64> { mle.iter().map(|x| c * x.abs()).collect() };
    println!("exact   {:.6?}", grad_exact(&mle, false));
    println!("abs_1e3 {:.6?}", fd_grad(loglik, &mle, &[1e-3; 3], true));
    println!("rel_1e3 {:.6?}", fd_grad(loglik, &mle, &scaled(1e-3), true));
    println!("rel_6e6 {:.6?}", fd_grad(loglik, &mle, &scaled(6e-6), true));
    // (c) a bug: dN/dr without the factor t
    let fd = fd_grad(loglik, &th, &hrel(6e-6), true);
    let gb = grad_exact(&th, true);
    let rb: Vec<String> = (0..3).map(|j| format!("{:.3e}", (gb[j] - fd[j]).abs() / fd[j].abs())).collect();
    println!("bug {:?}", rb);
}
