fn main() {
    // (a) sqrt(x + 1) - sqrt(x), naive vs rewritten
    for k in (2..=16).step_by(2) {
        let x = 10f64.powi(k);
        let naive = (x + 1.0).sqrt() - x.sqrt();
        let stable = 1.0 / ((x + 1.0).sqrt() + x.sqrt());
        // clamp to [0, 16]; `+ 0.0` turns a printed -0.0 into 0.0
        let dig = (-((naive - stable).abs() / stable).log10()).clamp(0.0, 16.0) + 0.0;
        println!("{:>4.0e}  {:>12.6e}  {:>12.6e}  {:5.1}", x, naive, stable, dig);
    }
    println!();

    // (b) small root of C^2 - S C + P L = 0
    let (p, l) = (0.01, 0.01);
    for k in 1..=7 {
        let kd = 10f64.powi(k);
        let s = p + l + kd;
        let c_naive = (s - (s * s - 4.0 * p * l).sqrt()) / 2.0;
        let c_stable = 2.0 * p * l / (s + (s * s - 4.0 * p * l).sqrt());
        println!("{:>3.0e}  {:>12.6e}  {:>12.6e}  {:>7.1e}", kd, c_naive, c_stable,
                 (c_naive - c_stable).abs() / c_stable);
    }
}
