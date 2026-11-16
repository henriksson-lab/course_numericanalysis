// Quiz: extrapolating the degree-4 interpolant beyond the last reading.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

fn lagrange(x: &[f64], y: &[f64], z: f64) -> f64 {
    (0..x.len())
        .map(|i| {
            let li: f64 = (0..x.len()).filter(|&j| j != i)
                .map(|j| (z - x[j]) / (x[i] - x[j])).product();
            li * y[i]
        })
        .sum()
}

fn main() {
    let t = [0.0, 2.0, 4.0, 6.0, 8.0];
    let y: Vec<f64> = t.iter().map(|&ti| od(ti)).collect();
    let p: Vec<f64> = [10.0, 12.0, 24.0].iter().map(|&z| lagrange(&t, &y, z)).collect();
    let truth: Vec<f64> = [10.0, 12.0, 24.0].iter().map(|&z| od(z)).collect();
    println!("{:?}", p);
    println!("{:?}", truth);
}
