// C2: a general Gillespie SSA; constitutive gene vs telegraph (bursty) gene.
use rand::{rngs::StdRng, RngExt, SeedableRng};

// state at the given output times; s[i][j] = change of species i by reaction j
fn ssa(x0: &[f64], s: &[Vec<f64>], propensity: &dyn Fn(&[f64]) -> Vec<f64>,
       times: &[f64], rng: &mut StdRng) -> Vec<Vec<f64>> {
    let mut x = x0.to_vec();
    let mut out = Vec::new();
    let (mut t, mut i) = (0.0, 0);
    loop {
        let a = propensity(&x);
        let a0: f64 = a.iter().sum();
        t += -(1.0 - rng.random::<f64>()).ln() / a0; // inf if a0 = 0
        while i < times.len() && times[i] < t {
            out.push(x.clone());
            i += 1;
        }
        if i == times.len() {
            return out;
        }
        let r = rng.random::<f64>() * a0; // reaction j with probability a_j / a0
        let mut j = 0;
        let mut cum = a[0];
        while cum <= r {
            j += 1;
            cum += a[j];
        }
        for (xi, si) in x.iter_mut().zip(s) {
            *xi += si[j];
        }
    }
}

fn mean_fano(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    let var = v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0);
    (m, var / m)
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    // (a) constitutive gene, species (m, p)
    let s2 = vec![vec![1.0, -1.0, 0.0, 0.0], vec![0.0, 0.0, 1.0, -1.0]];
    let gene = |x: &[f64]| vec![10.0, 1.0 * x[0], 2.0 * x[0], 0.2 * x[1]];
    let cells: Vec<Vec<f64>> = (0..1000)
        .map(|_| ssa(&[0.0, 0.0], &s2, &gene, &[30.0], &mut rng)[0].clone()).collect();
    let m: Vec<f64> = cells.iter().map(|c| c[0]).collect();
    let p: Vec<f64> = cells.iter().map(|c| c[1]).collect();
    let ((mm, fm), (mp, fp)) = (mean_fano(&m), mean_fano(&p));
    println!("constitutive: mean m {:.3}  mean p {:.3}  Fano m {:.4}  Fano p {:.4}", mm, mp, fm, fp);

    // (b) telegraph, species (g, m, p)
    let s3 = vec![vec![1.0, -1.0, 0.0, 0.0, 0.0, 0.0],
                  vec![0.0, 0.0, 1.0, -1.0, 0.0, 0.0],
                  vec![0.0, 0.0, 0.0, 0.0, 1.0, -1.0]];
    let telegraph = |x: &[f64]| vec![0.5 * (1.0 - x[0]), 4.5 * x[0], 100.0 * x[0],
                                     1.0 * x[1], 2.0 * x[1], 0.2 * x[2]];
    let cells: Vec<Vec<f64>> = (0..1000)
        .map(|_| ssa(&[0.0, 0.0, 0.0], &s3, &telegraph, &[30.0], &mut rng)[0].clone()).collect();
    let mb: Vec<f64> = cells.iter().map(|c| c[1]).collect();
    let pb: Vec<f64> = cells.iter().map(|c| c[2]).collect();
    let ((mm, fm), (mp, fp)) = (mean_fano(&mb), mean_fano(&pb));
    println!("bursty:       mean m {:.3}  mean p {:.3}  Fano m {:.4}  Fano p {:.4}", mm, mp, fm, fp);
    let zero = |v: &[f64]| v.iter().filter(|&&x| x == 0.0).count() as f64 / v.len() as f64;
    println!("P(m = 0): bursty {}  constitutive {}  Poisson(10) {:.4e}", zero(&mb), zero(&m),
             (-10f64).exp());
}
