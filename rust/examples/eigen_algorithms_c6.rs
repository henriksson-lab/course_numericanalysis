// C6: PCA of a simulated 20 000-cell x 2 000-gene sparse matrix, "full" route:
// form the 2000 x 2000 matrix Xc^T Xc from the sparse rows, then a dense
// symmetric eigendecomposition (nalgebra). A dense SVD of the 20000 x 2000 matrix
// in nalgebra (no LAPACK) is too slow to be worth waiting for.
// The random numbers differ from R and Python, so the numbers differ too.
use nalgebra::{DMatrix, SymmetricEigen};
use rand::{rngs::StdRng, RngExt, SeedableRng};
use rand_distr::{Distribution, Gamma, LogNormal, Normal, Poisson};
use sprs::CsMat;
use std::time::Instant;

fn simulate_cells(n_cells: usize, n_genes: usize, n_types: usize, seed: u64) -> (CsMat<f64>, Vec<usize>) {
    let mut rng = StdRng::seed_from_u64(seed);
    let gamma = Gamma::new(0.3, 0.5).unwrap();                  // mean count per gene
    let base: Vec<f64> = (0..n_genes).map(|_| gamma.sample(&mut rng)).collect();
    let normal = Normal::new(0.0, 1.0).unwrap();                // cell-type effects
    let fold: Vec<f64> = (0..n_genes * n_types).map(|_| 2f64.powf(normal.sample(&mut rng))).collect();
    let ctype: Vec<usize> = (0..n_cells).map(|_| rng.random_range(0..n_types)).collect();
    let lognorm = LogNormal::new(0.0, 0.3).unwrap();            // library size factor
    let lib: Vec<f64> = (0..n_cells).map(|_| lognorm.sample(&mut rng)).collect();
    // build the cells x genes matrix row by row in CSR form, storing only non-zeros
    let (mut indptr, mut indices, mut data) = (vec![0usize], Vec::new(), Vec::new());
    for i in 0..n_cells {
        for g in 0..n_genes {
            let lam = base[g] * fold[g * n_types + ctype[i]] * lib[i];
            let count: f64 = Poisson::new(lam).unwrap().sample(&mut rng);
            if count > 0.0 {
                indices.push(g);
                data.push(count.ln_1p());
            }
        }
        indptr.push(indices.len());
    }
    (CsMat::new((n_cells, n_genes), indptr, indices, data), ctype)
}

fn main() {
    let (x, _ctype) = simulate_cells(20000, 2000, 10, 1);
    let (n, p) = x.shape();
    let mut mu = vec![0.0; p];                                   // column means
    for row in x.outer_iterator() {
        for (j, &v) in row.iter() { mu[j] += v / n as f64; }
    }
    let mb = (x.nnz() * (8 + 8) + (n + 1) * 8) as f64 / 1e6;
    println!("{:.4} {:.1} {:.1}", x.nnz() as f64 / (n * p) as f64, mb, (8 * n * p) as f64 / 1e6);
    // covariance route: Xc^T Xc = X^T X - n mu mu^T, X^T X summed over sparse rows
    let t0 = Instant::now();
    let mut c = DMatrix::<f64>::zeros(p, p);
    for row in x.outer_iterator() {
        let nz: Vec<(usize, f64)> = row.iter().map(|(j, &v)| (j, v)).collect();
        for &(j, a) in &nz {
            for &(k, b) in &nz { c[(j, k)] += a * b; }
        }
    }
    for j in 0..p {
        for k in 0..p { c[(j, k)] -= n as f64 * mu[j] * mu[k]; }
    }
    let t_form = t0.elapsed().as_secs_f64();
    let t0 = Instant::now();
    let mut e: Vec<f64> = SymmetricEigen::new(c).eigenvalues.iter().cloned().collect();
    let t_eig = t0.elapsed().as_secs_f64();
    e.sort_by(|a, b| b.partial_cmp(a).unwrap());
    println!("form C {:.2} s   eigen {:.2} s", t_form, t_eig);
    let d: Vec<String> = e[..13].iter().map(|v| format!("{:.4}", v.sqrt())).collect();
    println!("{}", d.join(" "));
}
