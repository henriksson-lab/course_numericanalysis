// C5: long-run behaviour of the Ricker model for six growth rates.
const K: f64 = 100.0;

fn ricker(n: f64, r: f64) -> f64 { n * (r * (1.0 - n / K)).exp() }

fn main() {
    for r in [0.5, 1.5, 1.9, 2.3, 2.6, 3.0] {
        let mut n = 20.0;
        for _ in 0..1000 { n = ricker(n, r); } // throw away the transient
        let mut last = Vec::new();
        for _ in 0..64 { n = ricker(n, r); last.push(n); }
        // count distinct values after rounding to 6 decimals
        let mut keys: Vec<i64> = last.iter().map(|v| (v * 1e6).round() as i64).collect();
        keys.sort();
        keys.dedup();
        let first: Vec<String> = last[..4].iter().map(|v| format!("{}", (v * 100.0).round() / 100.0)).collect();
        println!("r = {r:.1}: g'(K) = {:4.1}, {:2} distinct values, e.g. {}",
                 1.0 - r, keys.len(), first.join(" "));
    }
    // No plotting here: print the cobweb corners (x_k, x_k), (x_k, x_k+1) and plot in R or Python.
}
