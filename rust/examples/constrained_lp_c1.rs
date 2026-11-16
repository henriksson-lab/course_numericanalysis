// C1: decay fit with background c, free and with c >= 0. Our crates have no
// L-BFGS-B, so this is hand-written: for fixed k the model A e^{-kt} + c is
// linear in (A, c), solved exactly; a golden-section search then finds k.
fn inner(k: f64, t: &[f64], y: &[f64], bounded: bool) -> (f64, f64, f64) {
    let e: Vec<f64> = t.iter().map(|&ti| (-k * ti).exp()).collect();
    let n = t.len() as f64;
    let (se, see) = (e.iter().sum::<f64>(), e.iter().map(|v| v * v).sum::<f64>());
    let sy: f64 = y.iter().sum();
    let sey: f64 = e.iter().zip(y).map(|(a, b)| a * b).sum();
    // normal equations [see se; se n] (A, c) = (sey, sy)
    let det = see * n - se * se;
    let (mut a, mut c) = ((sey * n - se * sy) / det, (see * sy - se * sey) / det);
    if bounded && c < 0.0 {
        c = 0.0; // the bound is active: refit A alone with c fixed at 0
        a = sey / see;
    }
    let sse = e.iter().zip(y).map(|(ei, yi)| (yi - a * ei - c).powi(2)).sum();
    (a, c, sse)
}

fn golden(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64) -> f64 {
    let g = (5f64.sqrt() - 1.0) / 2.0;
    while hi - lo > 1e-10 {
        let (x1, x2) = (hi - g * (hi - lo), lo + g * (hi - lo));
        if f(x1) < f(x2) { hi = x2 } else { lo = x1 }
    }
    (lo + hi) / 2.0
}

fn main() {
    let t: Vec<f64> = (0..9).map(|i| i as f64).collect();
    let y = [1.020, 0.577, 0.378, 0.243, 0.125, 0.062, 0.020, 0.010, -0.012];
    for (name, bounded) in [("free", false), ("bound", true)] {
        let k = golden(|k| inner(k, &t, &y, bounded).2, 0.05, 2.0);
        let (a, c, sse) = inner(k, &t, &y, bounded);
        println!("{:5} A {:.4}  k {:.4}  c {:.4}  SSE {:.4}  half-life {:.3} h",
                 name, a, k, c, sse, 2f64.ln() / k);
        if bounded {
            // dSSE/dc = -2 * sum of residuals
            let g: f64 = t.iter().zip(&y)
                .map(|(&ti, &yi)| -2.0 * (yi - a * (-k * ti).exp() - c)).sum();
            println!("dSSE/dc at the bounded fit: {:.4}", g);
        }
    }
}
