// C7: refit 500 bulk samples with fresh 5% noise: which estimates are stable?
// The random numbers differ from R and Python, so the last digits differ too.
use nalgebra::{DMatrix, DVector};
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, StandardNormal};

// The simulated signature matrix S (40 genes x 5 cell types), the true
// proportions, and the noise-free and noisy bulk profiles b0 and b.
fn signature() -> (DMatrix<f64>, DVector<f64>, DVector<f64>, DVector<f64>) {
    let mut s = DMatrix::from_element(40, 5, 10.0);
    for i in 0..16 { s[(i, 0)] = 100.0; s[(i, 1)] = 100.0; }  // T-cell genes
    for i in 16..40 { s[(i, 2 + (i - 16) / 8)] = 100.0; }     // B, NK, Mono markers
    for i in 0..40 {
        for j in 0..5 {                                        // gene-to-gene variation
            s[(i, j)] *= 1.0 + 0.5 * (((i + 1) * (j + 1)) as f64).sin();
        }
    }
    for i in 0..40 { s[(i, 1)] = s[(i, 0)]; }                  // CD8 = CD4 ...
    s[(0, 1)] = 0.98 * s[(0, 0)];                              // ... except gene 1
    s[(1, 1)] = 1.02 * s[(1, 0)];                              //     and gene 2
    let p_true = DVector::from_vec(vec![0.30, 0.15, 0.15, 0.10, 0.30]);
    let b0 = &s * &p_true;
    let b = DVector::from_fn(40, |i, _| b0[i] * (1.0 + 0.05 * (37.0 * (i + 1) as f64).sin()));
    (s, p_true, b0, b)
}

// least squares by Householder QR (as in the previous unit): R x = Q^T b
fn lstsq_qr(a: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let n = a.ncols();
    let qr = a.clone().qr();
    let mut qtb = b.clone();
    qr.q_tr_mul(&mut qtb);
    qr.r().solve_upper_triangular(&qtb.rows(0, n).into_owned()).unwrap()
}

// ridge (Tikhonov) through the SVD: filter factors f_i = sigma_i^2 / (sigma_i^2 + lambda^2)
fn ridge(a: &DMatrix<f64>, b: &DVector<f64>, lambda: f64) -> DVector<f64> {
    let svd = a.clone().svd(true, true);
    let u = svd.u.unwrap();
    let vt = svd.v_t.unwrap();
    let mut x = DVector::zeros(a.ncols());
    for (i, &s) in svd.singular_values.iter().enumerate() {
        let f = s * s / (s * s + lambda * lambda);
        x += f * u.column(i).dot(b) / s * vt.row(i).transpose();
    }
    x
}

// truncated SVD with the k largest singular values
fn tsvd(a: &DMatrix<f64>, b: &DVector<f64>, k: usize) -> DVector<f64> {
    let svd = a.clone().svd(true, true);
    let (u, vt) = (svd.u.unwrap(), svd.v_t.unwrap());
    let mut x = DVector::zeros(a.ncols());
    for i in 0..k {
        x += u.column(i).dot(b) / svd.singular_values[i] * vt.row(i).transpose();
    }
    x
}

// Non-negative least squares by the Lawson-Hanson active-set method (hand-written;
// R: nnls::nnls, Python: scipy.optimize.nnls). "free" = allowed to be positive.
fn nnls(a: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let n = a.ncols();
    let mut x = DVector::zeros(n);
    let mut free = vec![false; n];
    let tol = 1e-10 * (a.transpose() * b).abs().max();
    for _outer in 0..3 * n {
        let w = a.transpose() * (b - a * &x);           // -1/2 gradient of ||Ax - b||^2
        // the held variable whose increase reduces the residual fastest
        let j = (0..n).filter(|&j| !free[j] && w[j] > tol)
                      .max_by(|&i, &k| w[i].partial_cmp(&w[k]).unwrap());
        let Some(j) = j else { break };                 // none: optimal (KKT holds)
        free[j] = true;
        loop {
            // least squares on the free variables only
            let idx: Vec<usize> = (0..n).filter(|&i| free[i]).collect();
            let zf = lstsq_qr(&a.select_columns(&idx), b);
            let mut z = DVector::zeros(n);
            for (k, &i) in idx.iter().enumerate() { z[i] = zf[k]; }
            if idx.iter().all(|&i| z[i] > 0.0) { x = z; break; }
            // step from x towards z until the first free variable hits 0
            let alpha = idx.iter().filter(|&&i| z[i] <= 0.0)
                           .map(|&i| x[i] / (x[i] - z[i]))
                           .fold(f64::INFINITY, f64::min);
            x += alpha * (&z - &x);
            for &i in &idx {
                if x[i] <= 1e-14 { x[i] = 0.0; free[i] = false; }   // hold it at 0
            }
        }
    }
    x
}

fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    (m, (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0)).sqrt())
}

fn main() {
    let (s, _p_true, b0, _b) = signature();
    let mut rng = StdRng::seed_from_u64(1);
    // est[method][quantity] = 500 values; quantity 0 = CD4, 1 = CD8, 2 = CD4 + CD8
    let mut est = vec![vec![Vec::new(); 3]; 4];
    for _ in 0..500 {
        let br = DVector::from_fn(40, |i, _| {
            let z: f64 = StandardNormal.sample(&mut rng);
            b0[i] * (1.0 + 0.05 * z)                     // fresh 5% noise
        });
        let xs = [lstsq_qr(&s, &br), tsvd(&s, &br, 4), ridge(&s, &br, 10.0), nnls(&s, &br)];
        for (m, x) in xs.iter().enumerate() {
            est[m][0].push(x[0]);
            est[m][1].push(x[1]);
            est[m][2].push(x[0] + x[1]);
        }
    }
    println!("      ls  tsvd ridge  nnls");
    for (q, name) in ["CD4", "CD8", "sum"].iter().enumerate() {
        let ms: Vec<(f64, f64)> = (0..4).map(|m| mean_sd(&est[m][q])).collect();
        let means: Vec<String> = ms.iter().map(|p| format!("{:5.3}", p.0)).collect();
        let sds: Vec<String> = ms.iter().map(|p| format!("{:5.3}", p.1)).collect();
        println!("{} {}  | sd {}", name, means.join(" "), sds.join(" "));
    }
    let zero = |q: usize| est[3][q].iter().filter(|&&x| x == 0.0).count() as f64 / 500.0;
    println!("{} {}", zero(1), zero(0));
}
