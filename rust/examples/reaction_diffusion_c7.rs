// C7: a 2D Turing pattern (Schnakenberg), IMEX Euler with FFT for diffusion.
use rustfft::{num_complex::Complex, FftPlanner};

/// Park-Miller generator, the same numbers as in R and Python.
fn lcg(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..n).map(|_| { s = (16807 * s) % 2147483647; s as f64 / 2147483647.0 }).collect()
}

/// 2D FFT of an N x N grid stored column by column (w[i + n*j]): rows, then columns.
fn fft2(w: &mut [Complex<f64>], n: usize, planner: &mut FftPlanner<f64>, inverse: bool) {
    let fft = if inverse { planner.plan_fft_inverse(n) } else { planner.plan_fft_forward(n) };
    for col in w.chunks_mut(n) { fft.process(col); } // along i, for each j
    let mut tmp = vec![Complex::new(0.0, 0.0); n];
    for i in 0..n { // along j, for each i
        for j in 0..n { tmp[j] = w[i + n * j]; }
        fft.process(&mut tmp);
        for j in 0..n { w[i + n * j] = tmp[j]; }
    }
}

/// Schnakenberg on an N x N periodic grid of side 80: IMEX Euler, the implicit
/// five-point diffusion solved by FFT. Returns u (column by column).
fn schnak(n: usize, dt: f64, tend: f64, d: f64) -> Vec<f64> {
    let (a, b, lbox) = (0.1, 0.9, 80.0);
    let dx = lbox / n as f64;
    let s: Vec<f64> = (0..n)
        .map(|i| 4.0 / (dx * dx) * (std::f64::consts::PI * i as f64 / n as f64).sin().powi(2))
        .collect();
    let kap: Vec<f64> = (0..n * n).map(|idx| s[idx % n] + s[idx / n]).collect();
    let z = lcg(2 * n * n, 1);
    let mut u: Vec<f64> = (0..n * n).map(|i| 1.0 + 0.02 * (z[i] - 0.5)).collect();
    let mut v: Vec<f64> = (0..n * n).map(|i| 0.9 + 0.02 * (z[n * n + i] - 0.5)).collect();
    let mut planner = FftPlanner::new();
    let mut solve_diff = |w: &[f64], dw: f64| -> Vec<f64> {
        let mut c: Vec<Complex<f64>> = w.iter().map(|&x| Complex::new(x, 0.0)).collect();
        fft2(&mut c, n, &mut planner, false);
        for (ci, k) in c.iter_mut().zip(&kap) { *ci /= 1.0 + dt * dw * k; }
        fft2(&mut c, n, &mut planner, true);
        c.iter().map(|x| x.re / (n * n) as f64).collect()
    };
    for _ in 0..(tend / dt).round() as usize {
        let mut ru = vec![0.0; n * n];
        let mut rv = vec![0.0; n * n];
        for i in 0..n * n {
            let uuv = u[i] * u[i] * v[i];
            ru[i] = u[i] + dt * (a - u[i] + uuv); // reaction: explicit
            rv[i] = v[i] + dt * (b - uuv);
        }
        u = solve_diff(&ru, 1.0); // diffusion: implicit
        v = solve_diff(&rv, d);
    }
    u
}

/// Number of local maxima (8 neighbours, periodic) above the mean.
fn spots(u: &[f64], n: usize) -> usize {
    let mean = u.iter().sum::<f64>() / (n * n) as f64;
    let mut count = 0;
    for j in 0..n {
        for i in 0..n {
            let c = u[i + n * j];
            let mut ismax = c > mean;
            for di in [n - 1, 0, 1] {
                for dj in [n - 1, 0, 1] {
                    if di == 0 && dj == 0 { continue; }
                    if !(c > u[(i + di) % n + n * ((j + dj) % n)]) { ismax = false; }
                }
            }
            if ismax { count += 1; }
        }
    }
    count
}

fn main() {
    let t0 = std::time::Instant::now();
    let n = 80;
    let u = schnak(n, 0.1, 200.0, 20.0);
    println!("{:.2} s", t0.elapsed().as_secs_f64());
    let (umin, umax) = u.iter().fold((f64::MAX, f64::MIN), |m, &x| (m.0.min(x), m.1.max(x)));
    let mean = u.iter().sum::<f64>() / (n * n) as f64;
    println!("spots {}, u from {:.6} to {:.6}, mean {:.6}", spots(&u, n), umin, umax, mean);
}
