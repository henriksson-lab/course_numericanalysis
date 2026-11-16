// C4: growth-curve data from the slides, the logistic model and the SSE.
fn logistic(p: &[f64], t: f64) -> f64 {
    let (r, k, n0) = (p[0], p[1], p[2]); // p = (r, K, N0)
    k / (1.0 + (k / n0 - 1.0) * (-r * t).exp())
}

fn main() {
    let t_min: Vec<f64> = (0..=12).map(|i| 120.0 * i as f64).collect(); // minutes
    let od = [0.001, 0.047, 0.139, 0.282, 0.581, 0.859, 1.049,
              1.161, 1.152, 1.217, 1.182, 1.176, 1.185];
    // A closure: a small function that "captures" t_min and od from above
    let sse = |p: &[f64]| -> f64 {
        t_min.iter().zip(&od).map(|(&t, &y)| (y - logistic(p, t)).powi(2)).sum()
    };
    println!("SSE at (0.5, 1, 0.05): {:.4}", sse(&[0.5, 1.0, 0.05]));
}
