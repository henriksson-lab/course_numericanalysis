// C7: matrix-free Lanczos (full reorthogonalisation) for the top eigenpairs of
// Xc^T Xc, with X sparse (sprs) and the centring done implicitly. Rust has no
// counterpart to irlba / eigsh in our crate list, so Lanczos is hand-written.
// The random numbers differ from R and Python, so the numbers differ too.
use nalgebra::{DMatrix, DVector, SymmetricEigen};
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

// Xc^T Xc v without forming Xc:  w = X v - 1 (mu . v),  out = X^T w - mu (1 . w)
fn cv(x: &CsMat<f64>, mu: &DVector<f64>, v: &DVector<f64>) -> DVector<f64> {
    let mv = mu.dot(v);
    let w: Vec<f64> = x.outer_iterator()
        .map(|row| row.iter().map(|(j, &a)| a * v[j]).sum::<f64>() - mv)
        .collect();
    let wsum: f64 = w.iter().sum();
    let mut out = DVector::<f64>::zeros(v.len());
    for (i, row) in x.outer_iterator().enumerate() {
        for (j, &a) in row.iter() { out[j] += a * w[i]; }
    }
    out - mu * wsum
}

// Lanczos with full reorthogonalisation: Ritz values (decreasing) and vectors
fn lanczos_eig(afun: &dyn Fn(&DVector<f64>) -> DVector<f64>, n: usize, m: usize) -> (Vec<f64>, DMatrix<f64>) {
    let mut qm = DMatrix::<f64>::zeros(n, m);
    let (mut alpha, mut beta) = (vec![0.0; m], vec![0.0; m]);
    let mut q = DVector::from_element(n, 1.0 / (n as f64).sqrt());
    let mut q_old = DVector::<f64>::zeros(n);
    let mut b = 0.0;
    for j in 0..m {
        qm.set_column(j, &q);
        let mut w = afun(&q) - &q_old * b;
        alpha[j] = w.dot(&q);
        w -= &q * alpha[j];
        let qj = qm.columns(0, j + 1);
        w -= &qj * (qj.transpose() * &w);
        b = w.norm();
        beta[j] = b;
        q_old = q;
        q = w / b;
    }
    let t = DMatrix::from_fn(m, m, |i, j| {
        if i == j { alpha[i] } else if i == j + 1 || j == i + 1 { beta[i.min(j)] } else { 0.0 }
    });
    let eig = SymmetricEigen::new(t);
    let mut idx: Vec<usize> = (0..m).collect();
    idx.sort_by(|&a, &b| eig.eigenvalues[b].partial_cmp(&eig.eigenvalues[a]).unwrap());
    let theta: Vec<f64> = idx.iter().map(|&i| eig.eigenvalues[i]).collect();
    let s = DMatrix::from_fn(m, m, |r, c| eig.eigenvectors[(r, idx[c])]);
    (theta, qm * s)
}

fn main() {
    let (x, _ctype) = simulate_cells(20000, 2000, 10, 1);
    let (n, p) = x.shape();
    let mut mu = DVector::<f64>::zeros(p);
    for row in x.outer_iterator() {
        for (j, &v) in row.iter() { mu[j] += v / n as f64; }
    }
    let afun = |v: &DVector<f64>| cv(&x, &mu, v);
    let t0 = Instant::now();
    let (theta, y) = lanczos_eig(&afun, p, 40);
    let t_l = t0.elapsed().as_secs_f64();
    let d: Vec<String> = theta[..12].iter().map(|v| format!("{:.4}", v.sqrt())).collect();
    println!("{}", d.join(" "));
    let res: Vec<String> = (0..12).map(|i| {               // the check: ||C y - theta y|| / theta
        let yi = y.column(i).into_owned();
        format!("{:.1e}", (afun(&yi) - &yi * theta[i]).norm() / theta[i])
    }).collect();
    println!("{}", res.join(" "));
    println!("lanczos {:.2} s", t_l);
}
