// C2: measure the order from the slope, and with the ratio test.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

// Least-squares line y = a + b x; returns (intercept, slope), like coef(lm(y ~ x)).
fn fit_line(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let sxy: f64 = x.iter().zip(y).map(|(xi, yi)| (xi - mx) * (yi - my)).sum();
    let sxx: f64 = x.iter().map(|xi| (xi - mx).powi(2)).sum();
    let b = sxy / sxx;
    (my - b * mx, b)
}

fn main() {
    let t0 = 3.0;
    let exact = R * od(t0) * (1.0 - od(t0) / K);

    let h: Vec<f64> = (0..7).map(|k| 10f64.powf(-(1.0 + 0.5 * k as f64))).collect();
    let fwd: Vec<f64> = h.iter().map(|&h| ((od(t0 + h) - od(t0)) / h - exact).abs()).collect();
    let cen: Vec<f64> = h.iter().map(|&h| ((od(t0 + h) - od(t0 - h)) / (2.0 * h) - exact).abs()).collect();
    let log10 = |v: &[f64]| v.iter().map(|x| x.log10()).collect::<Vec<f64>>();
    println!("{:?}", fit_line(&log10(&h), &log10(&fwd))); // (intercept, slope)
    println!("{:?}", fit_line(&log10(&h), &log10(&cen)));

    let d = |h: f64| (od(t0 + h) - od(t0 - h)) / (2.0 * h); // a closure: a small inline function
    let ratios = |h0: f64| -> Vec<f64> {
        let dv: Vec<f64> = (0..6).map(|k| d(h0 / 2f64.powi(k))).collect();
        let diff: Vec<f64> = dv.windows(2).map(|w| w[0] - w[1]).collect(); // D(h) - D(h/2), ...
        diff.windows(2).map(|w| w[0] / w[1]).collect() // should be 2^p
    };
    println!("{:?}", ratios(0.1));
    println!("{:?}", ratios(1e-8));
}
