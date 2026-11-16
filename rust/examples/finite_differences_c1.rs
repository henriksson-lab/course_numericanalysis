// C1: four stencils, their errors, best step and order.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

// slope of the least-squares line through (x, y)
fn slope(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    sxy / sxx
}

fn main() {
    let t0 = 3.0;
    let d1 = R * od(t0) * (1.0 - od(t0) / K);
    let d2 = R * d1 * (1.0 - 2.0 * od(t0) / K);
    let h: Vec<f64> = (1..=12).map(|k| 10f64.powi(-k)).collect();
    let names = ["fwd", "cen", "one", "sec"];
    let err: Vec<Vec<f64>> = vec![
        h.iter().map(|&h| ((od(t0 + h) - od(t0)) / h - d1).abs()).collect(),
        h.iter().map(|&h| ((od(t0 + h) - od(t0 - h)) / (2.0 * h) - d1).abs()).collect(),
        h.iter().map(|&h| ((-3.0 * od(t0) + 4.0 * od(t0 + h) - od(t0 + 2.0 * h)) / (2.0 * h) - d1).abs()).collect(),
        h.iter().map(|&h| ((od(t0 + h) - 2.0 * od(t0) + od(t0 - h)) / (h * h) - d2).abs()).collect(),
    ];
    println!("    h      fwd      cen      one      sec");
    for i in 0..h.len() {
        print!("{:.0e}", h[i]);
        for e in &err { print!("  {:.2e}", e[i]); }
        println!();
    }
    print!("best h:");
    for (k, e) in err.iter().enumerate() {
        let i = (0..h.len()).min_by(|&a, &b| e[a].partial_cmp(&e[b]).unwrap()).unwrap();
        print!(" {} {:.0e}", names[k], h[i]);
    }
    println!();
    let lh: Vec<f64> = h[..3].iter().map(|x| x.log10()).collect(); // h >= 1e-3
    print!("slope: ");
    for (k, e) in err.iter().enumerate() {
        let le: Vec<f64> = e[..3].iter().map(|x| x.log10()).collect();
        print!(" {} {:.4}", names[k], slope(&lh, &le));
    }
    println!();
}
