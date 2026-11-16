// C1: a toy LCG, the RANDU planes, and reproducible seeds.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use std::collections::HashSet;

fn lcg(n: usize, seed: u64, a: u64, c: u64, m: u64) -> Vec<u64> {
    let mut x = vec![seed];
    for _ in 1..n {
        let last = *x.last().unwrap();
        x.push((a * last + c) % m);
    }
    x
}

fn n_distinct(v: &[f64]) -> usize {
    // distinct values after rounding to 3 decimals
    let set: HashSet<i64> = v.iter().map(|x| (x * 1000.0).round() as i64).collect();
    set.len()
}

fn main() {
    println!("{:?}", lcg(17, 1, 5, 1, 16));
    println!("{:?}", lcg(6, 1, 4, 1, 16));

    // RANDU: x_{n+1} = 65539 x_n mod 2^31 (u64 arithmetic is exact here)
    let u: Vec<f64> = lcg(30000, 1, 65539, 0, 1 << 31)
        .iter().map(|&x| x as f64 / 2f64.powi(31)).collect();
    let v: Vec<f64> = (0..29998).map(|i| 9.0 * u[i] - 6.0 * u[i + 1] + u[i + 2]).collect();
    println!("{}", n_distinct(&v));
    let mut rng = StdRng::seed_from_u64(1);
    let w: Vec<f64> = (0..30000).map(|_| rng.random::<f64>()).collect();
    let v: Vec<f64> = (0..29998).map(|i| 9.0 * w[i] - 6.0 * w[i + 1] + w[i + 2]).collect();
    println!("{}", n_distinct(&v));

    let mut r1 = StdRng::seed_from_u64(42);
    let mut r2 = StdRng::seed_from_u64(42);
    let a: Vec<f64> = (0..3).map(|_| r1.random::<f64>()).collect();
    let b: Vec<f64> = (0..3).map(|_| r2.random::<f64>()).collect();
    println!("{}", a == b);
}
