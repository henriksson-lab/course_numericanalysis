// C4: doubling experiment for the gene-gene correlation matrix.
// Run in release mode, or the timings mean nothing:
//   cargo run -q --release --example convergence_cost_c4
// The correlation is a hand-written loop over all gene pairs (no BLAS, one core),
// so the n^2 scaling shows cleanly. Times depend on your machine; the
// random counts differ from R and Python.
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, Poisson};
use std::hint::black_box;
use std::time::Instant;

// Correlation of the columns of a cells x genes matrix, stored as one Vec per gene.
fn cor(genes: &[Vec<f64>]) -> Vec<Vec<f64>> {
    // centre each gene and scale it to length 1; then correlation = dot product
    let z: Vec<Vec<f64>> = genes.iter().map(|g| {
        let m = g.iter().sum::<f64>() / g.len() as f64;
        let c: Vec<f64> = g.iter().map(|v| v - m).collect();
        let norm = c.iter().map(|v| v * v).sum::<f64>().sqrt();
        c.iter().map(|v| v / norm).collect()
    }).collect();
    let n = z.len();
    let mut r = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in i..n { // each pair once: n(n+1)/2 dot products of length `cells`
            let s: f64 = z[i].iter().zip(&z[j]).map(|(a, b)| a * b).sum();
            r[i][j] = s;
            r[j][i] = s;
        }
    }
    r
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let pois = Poisson::new(2.0).unwrap();
    let cells = 1000;
    let genes = [500, 1000, 2000, 4000];
    let mut secs = vec![0.0; genes.len()];
    for (i, &g) in genes.iter().enumerate() {
        let x: Vec<Vec<f64>> = (0..g)
            .map(|_| (0..cells).map(|_| pois.sample(&mut rng)).collect())
            .collect(); // cells x genes, one Vec per gene
        let start = Instant::now();
        black_box(cor(black_box(&x)));
        secs[i] = start.elapsed().as_secs_f64();
    }
    println!("{:?}", secs.iter().map(|s| (s * 1000.0).round() / 1000.0).collect::<Vec<_>>());
    println!("{:?}", secs.windows(2).map(|w| (w[1] / w[0] * 100.0).round() / 100.0).collect::<Vec<_>>());
    // slope of log10(secs) against log10(genes), by least squares
    let lx: Vec<f64> = genes.iter().map(|&g| (g as f64).log10()).collect();
    let ly: Vec<f64> = secs.iter().map(|s| s.log10()).collect();
    let (mx, my) = (lx.iter().sum::<f64>() / 4.0, ly.iter().sum::<f64>() / 4.0);
    let slope = lx.iter().zip(&ly).map(|(a, b)| (a - mx) * (b - my)).sum::<f64>()
        / lx.iter().map(|a| (a - mx).powi(2)).sum::<f64>();
    println!("{slope}");
    println!("{}", secs[3] * (20000.0 / 4000.0_f64).powi(2) * (100000.0 / cells as f64) / 3600.0); // hours
}
