// C6: non-negative least squares, hand-written (Lawson-Hanson), and the optimality check.
use nalgebra::{DMatrix, DVector};

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

fn fmt(v: &DVector<f64>, d: usize) -> String {
    v.iter().map(|x| format!("{:.*}", d, x)).collect::<Vec<_>>().join(" ")
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

fn main() {
    // B5: rows (1, 1), (1, 1), (1, 1.1)
    let a = DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 1.0, 1.0, 1.0, 1.1]);
    let b5 = DVector::from_vec(vec![1.0, 1.0, 0.9]);
    println!("{}", fmt(&nnls(&a, &b5), 7));

    let (s, _p_true, _b0, b) = signature();
    let x = nnls(&s, &b);
    println!("{}", fmt(&x, 4));
    let w = s.transpose() * (&b - &s * &x);          // w = S^T (b - S x)
    println!("{}", fmt(&w, 6));
    let free_max = (0..5).filter(|&j| x[j] > 0.0).map(|j| w[j].abs()).fold(0.0, f64::max);
    println!("{:e}", free_max);                      // ~ 0 for the free variables
    println!("{}", (0..5).filter(|&j| x[j] == 0.0).all(|j| w[j] <= 1e-8));
}
