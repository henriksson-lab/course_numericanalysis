// Background: the PK data, the Bateman function and the residual sum of squares
const D: f64 = 100.0; // oral dose, mg
const TT: [f64; 12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0]; // h
const Y: [f64; 12] = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61]; // mg/L

fn bateman(t: f64, ka: f64, ke: f64, v: f64) -> f64 {
    D * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp())
}

fn rss(p: [f64; 3]) -> f64 {
    TT.iter().zip(Y.iter()).map(|(&t, &y)| (y - bateman(t, p[0], p[1], p[2])).powi(2)).sum()
}

fn main() {
    println!("{}", rss([1.0, 0.2, 20.0]));                       // at the simulation parameters
    // no nls in Rust: we write Gauss-Newton ourselves in C2; here the RSS at nls's estimate
    println!("{}", rss([1.0150721, 0.1906157, 20.2162104]));
}
