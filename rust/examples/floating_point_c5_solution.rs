use statrs::distribution::{Discrete, Poisson};

// index of the FIRST largest value, like R's which.max and np.argmax
fn which_max(v: &[f64]) -> usize {
    let mut best = 0;
    for i in 1..v.len() {
        if v[i] > v[best] {
            best = i;
        }
    }
    best
}

fn main() {
    let base: [u64; 10] = [2, 4, 3, 1, 5, 3, 2, 6, 3, 2];
    let counts: Vec<u64> = base.repeat(1000);                  // 10 000 cells
    let lambda: Vec<f64> = (0..=200).map(|i| 2.0 + i as f64 * 0.01).collect();

    // likelihood (product) and log-likelihood (sum) for one lambda
    let lik = |l: f64, x: &[u64]| {
        let p = Poisson::new(l).unwrap();
        x.iter().map(|&k| p.pmf(k)).product::<f64>()
    };
    let loglik = |l: f64, x: &[u64]| {
        let p = Poisson::new(l).unwrap();
        x.iter().map(|&k| p.ln_pmf(k)).sum::<f64>()
    };
    let liks: Vec<f64> = lambda.iter().map(|&l| lik(l, &counts)).collect();
    let logliks: Vec<f64> = lambda.iter().map(|&l| loglik(l, &counts)).collect();

    let min = liks.iter().cloned().fold(f64::INFINITY, f64::min);
    let max = liks.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("{} {}", min, max);
    println!("{:.2}", lambda[which_max(&liks)]);
    println!("{:.2}", lambda[which_max(&logliks)]);
    println!("{}", counts.iter().sum::<u64>() as f64 / counts.len() as f64);
    println!("{:e} {:e}", lik(3.0, &counts[..300]), lik(3.0, &counts[..400]));
}
