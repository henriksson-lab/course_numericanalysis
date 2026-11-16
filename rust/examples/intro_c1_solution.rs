// C1: errors of Taylor polynomials for e^x, their order, and a log-log slope.
fn main() {
    let x: [f64; 4] = [0.4, 0.2, 0.1, 0.05]; // a fixed-size array of 4 f64
    // err[i] holds the errors of P1, P2, P3 at x[i] (an array of 3 f64 per row)
    let err: Vec<[f64; 3]> = x
        .iter()
        .map(|&xi| {
            let p1 = 1.0 + xi;
            let p2 = p1 + xi * xi / 2.0;
            let p3 = p2 + xi.powi(3) / 6.0;
            [xi.exp() - p1, xi.exp() - p2, xi.exp() - p3]
        })
        .collect(); // collect the results into a Vec (a growable vector)
    for (xi, e) in x.iter().zip(&err) {
        println!("x = {:4.2}   {:.2e}  {:.2e}  {:.2e}", xi, e[0], e[1], e[2]);
    }

    // log2 of (error at x) / (error at x/2), one row per pair of consecutive x
    println!("        P1    P2    P3");
    for i in 0..3 {
        let r: Vec<f64> = (0..3).map(|j| (err[i][j] / err[i + 1][j]).log2()).collect();
        println!("[{}]   {:.2}  {:.2}  {:.2}", i + 1, r[0], r[1], r[2]);
    }

    // 20 points from 10^-2 to 10^0, evenly spaced on a log scale
    let xs: Vec<f64> = (0..20).map(|i| 10f64.powf(-2.0 + 2.0 * i as f64 / 19.0)).collect();
    let lx: Vec<f64> = xs.iter().map(|v| v.log10()).collect();
    let ly: Vec<f64> = xs.iter().map(|v| (v.exp() - (1.0 + v)).abs().log10()).collect();
    // No plotting here (R: plot, Python: matplotlib); the straight-line fit
    // that lm / polyfit do is written out by hand: slope = cov(x, y) / var(x).
    let n = lx.len() as f64;
    let (mx, my) = (lx.iter().sum::<f64>() / n, ly.iter().sum::<f64>() / n);
    let sxy: f64 = lx.iter().zip(&ly).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = lx.iter().map(|a| (a - mx).powi(2)).sum();
    let slope = sxy / sxx;
    println!("intercept {:.7}   slope {:.7}", my - slope * mx, slope);
}
