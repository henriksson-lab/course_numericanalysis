// C6: Fokker-Planck PDE (flux form, explicit Euler) against a histogram's moments of EM paths.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::{Distribution, Normal, StandardNormal};

const THETA: f64 = 1.0; // per ms
const D: f64 = 430.0; // nm^2/ms

// mass, mean and variance of the PDE solution at t_end
fn fp(t_end: f64, dx: f64, dt: f64, upwind: bool) -> [f64; 3] {
    let n = ((260.0 / dx).round() as usize) + 1;
    let x: Vec<f64> = (0..n).map(|i| -100.0 + i as f64 * dx).collect();
    let u: Vec<f64> = (0..n - 1).map(|i| -THETA * (x[i + 1] - dx / 2.0)).collect(); // at the faces
    let mut p: Vec<f64> = x.iter()
        .map(|&xi| (-(xi - 60.0).powi(2) / 50.0).exp() / (2.0 * std::f64::consts::PI * 25.0).sqrt())
        .collect();
    for _ in 0..(t_end / dt).round() as usize {
        let f: Vec<f64> = (0..n - 1).map(|i| {
            let pf = if upwind { if u[i] > 0.0 { p[i] } else { p[i + 1] } } else { (p[i] + p[i + 1]) / 2.0 };
            u[i] * pf - D * (p[i + 1] - p[i]) / dx // drift + Fick
        }).collect();
        for i in 0..n {
            let fr = if i < n - 1 { f[i] } else { 0.0 }; // no flux through the ends
            let fl = if i > 0 { f[i - 1] } else { 0.0 };
            p[i] -= dt / dx * (fr - fl);
        }
    }
    let mass: f64 = p.iter().sum::<f64>() * dx;
    let mean: f64 = x.iter().zip(&p).map(|(a, b)| a * b).sum::<f64>() * dx;
    let var: f64 = x.iter().zip(&p).map(|(a, b)| (a - mean).powi(2) * b).sum::<f64>() * dx;
    [mass, mean, var]
}

fn main() {
    let sigma = (2.0 * D).sqrt();
    let mut rng = StdRng::seed_from_u64(1);
    let start = Normal::new(60.0, 5.0).unwrap();
    let mut x: Vec<f64> = (0..10000).map(|_| start.sample(&mut rng)).collect();
    let h = 0.01;
    for i in 1..=300 {
        for xi in x.iter_mut() {
            let z: f64 = rng.sample(StandardNormal);
            *xi += -THETA * *xi * h + sigma * h.sqrt() * z;
        }
        if i == 100 || i == 300 {
            let t = i as f64 * h;
            let n = x.len() as f64;
            let m = x.iter().sum::<f64>() / n;
            let v = x.iter().map(|a| (a - m).powi(2)).sum::<f64>() / (n - 1.0);
            let ex = [1.0, 60.0 * (-THETA * t).exp(),
                      D / THETA * (1.0 - (-2.0 * THETA * t).exp()) + 25.0 * (-2.0 * THETA * t).exp()];
            println!("t = {}: mass, mean, var", t);
            println!("  pde     {:?}", fp(t, 2.0, 0.002, false).map(|v| (v * 1e4).round() / 1e4));
            println!("  upwind  {:?}", fp(t, 2.0, 0.002, true).map(|v| (v * 1e4).round() / 1e4));
            println!("  em      {:?}", [1.0, (m * 1e4).round() / 1e4, (v * 1e4).round() / 1e4]);
            println!("  exact   {:?}", ex.map(|v| (v * 1e4).round() / 1e4));
        }
    }
}
