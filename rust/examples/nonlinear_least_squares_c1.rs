// C1: the PK data and the Bateman function; tmax, Cmax, AUC; the 0/0 at ka = ke
const D: f64 = 100.0;
const TT: [f64; 12] = [0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 10.0, 12.0];
const Y: [f64; 12] = [1.02, 1.88, 2.34, 2.98, 3.27, 3.26, 3.17, 2.77, 1.92, 1.23, 1.00, 0.61];

fn bateman(t: f64, ka: f64, ke: f64, v: f64) -> f64 {
    D * ka / (v * (ka - ke)) * ((-ke * t).exp() - (-ka * t).exp())
}

fn rss(p: [f64; 3]) -> f64 {
    TT.iter().zip(Y.iter()).map(|(&t, &y)| (y - bateman(t, p[0], p[1], p[2])).powi(2)).sum()
}

fn main() {
    println!("{} {}", rss([1.0, 0.2, 20.0]), rss([0.5, 0.1, 10.0]));
    let tmax = (1.0f64 / 0.2).ln() / (1.0 - 0.2);
    println!("{} {} {}", tmax, bateman(tmax, 1.0, 0.2, 20.0), D / (20.0 * 0.2));
    for t in [1.0, 4.0] {
        println!("{} {} {}",
                 bateman(t, 0.2, 0.2, 20.0),                    // ka = ke: 0/0
                 D * 0.2 * t * (-0.2 * t).exp() / 20.0,         // the limit ka -> ke
                 bateman(t, 0.2 + 1e-9, 0.2, 20.0));            // close to ka = ke
    }
    // a table instead of a plot
    for t in [0.0, 1.0, 2.0, 4.0, 8.0, 12.0] {
        println!("t = {:>4}  C = {:.3}", t, bateman(t, 1.0, 0.2, 20.0));
    }
}
