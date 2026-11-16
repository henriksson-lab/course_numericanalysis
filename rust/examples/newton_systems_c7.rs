// C7: how many steady states for each promoter strength alpha?
use nalgebra::{DMatrix, DVector};

const N: f64 = 2.0;

fn f(x: &DVector<f64>, a: f64) -> DVector<f64> {
    let (u, v) = (x[0], x[1]);
    DVector::from_vec(vec![a / (1.0 + v.powf(N)) - u, a / (1.0 + u.powf(N)) - v])
}

fn jac(x: &DVector<f64>, a: f64) -> DMatrix<f64> {
    let (u, v) = (x[0], x[1]);
    DMatrix::from_row_slice(2, 2, &[
        -1.0, -a * N * v.powf(N - 1.0) / (1.0 + v.powf(N)).powi(2),
        -a * N * u.powf(N - 1.0) / (1.0 + u.powf(N)).powi(2), -1.0])
}

// last iterate and number of steps; None if J is exactly singular
fn newton_last(a: f64, x0: &DVector<f64>, tol: f64, maxit: usize) -> Option<(DVector<f64>, usize)> {
    let mut x = x0.clone();
    for k in 1..=maxit {
        let delta = jac(&x, a).lu().solve(&(-f(&x, a)))?;
        x += &delta;
        if delta.norm() < tol { return Some((x, k)); }
    }
    Some((x, maxit))
}

fn main() {
    let g: Vec<f64> = (0..8).map(|i| 0.25 + i as f64 * 4.75 / 7.0).collect();
    for i in 0..9 {
        let a = 1.0 + 0.5 * i as f64;
        let mut distinct: Vec<[f64; 2]> = Vec::new();
        let mut hit_maxit = 0;
        for &v in &g {
            for &u in &g {
                let Some((x, k)) = newton_last(a, &DVector::from_vec(vec![u, v]), 1e-12, 50) else { continue };
                if k == 50 { hit_maxit += 1; }
                if x.iter().all(|c| c.is_finite()) && f(&x, a).norm() < 1e-10 {
                    let key = [(x[0] * 1e6).round() / 1e6, (x[1] * 1e6).round() / 1e6];
                    if !distinct.contains(&key) { distinct.push(key); }
                }
            }
        }
        let umax = distinct.iter().map(|k| k[0]).fold(f64::MIN, f64::max);
        println!("alpha = {:.1}: {:2} distinct roots; largest u = {:.6}; {:2} runs hit maxit",
                 a, distinct.len(), umax, hit_maxit);
    }
    println!("{}", (5.0 + (5.0f64 * 5.0 - 4.0).sqrt()) / 2.0); // B5 formula at alpha = 5
}
