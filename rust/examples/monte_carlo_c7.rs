// C7: motif hits in promoters, simulated (permutation-style) vs exact binomial.
use rand::{rngs::StdRng, RngExt, SeedableRng};
use statrs::distribution::{Binomial, DiscreteCDF};

fn main() {
    let (g, l) = (1_000_000u64, 10u64);
    let pstart: Vec<u64> = (0..50).map(|j| 10001 + 20000 * j).collect();
    let pend: Vec<u64> = pstart.iter().map(|s| s + 999).collect();
    let n_hits = 200u64;
    let p = 50.0 * (1000 + l - 1) as f64 / (g - l + 1) as f64;

    // does a hit starting at s overlap a promoter?
    let overlaps = |s: u64| -> bool {
        let i = pstart.partition_point(|&ps| ps <= s + l - 1); // like findInterval
        i > 0 && pend[i - 1] >= s
    };

    let mut rng = StdRng::seed_from_u64(1);
    let b = 10000;
    let k: Vec<u64> = (0..b).map(|_| {
        (0..n_hits).filter(|_| overlaps(rng.random_range(1..=g - l + 1))).count() as u64
    }).collect();
    let binom = Binomial::new(p, n_hits).unwrap();
    for k_obs in [18u64, 25] {
        let c = k.iter().filter(|&&x| x >= k_obs).count();
        let ph = c as f64 / b as f64;
        println!("k_obs {}: simulated {} of {}, p = {:.5} (se {:.5}), (k+1)/(B+1) = {:.5}, exact {:.3e}",
                 k_obs, c, b, ph, (ph * (1.0 - ph) / b as f64).sqrt(),
                 (c + 1) as f64 / (b + 1) as f64, binom.sf(k_obs - 1)); // sf(x) = P(K > x)
    }
    let m = k.iter().sum::<u64>() as f64 / b as f64;
    let v = k.iter().map(|&x| (x as f64 - m).powi(2)).sum::<f64>() / (b - 1) as f64;
    println!("{} {} {} {}", m, n_hits as f64 * p, v, n_hits as f64 * p * (1.0 - p));
}
