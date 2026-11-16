// C4: least squares through the SVD, the Picard table, and the truncated SVD.
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

// truncated SVD: keep the k largest singular values
fn tsvd(a: &DMatrix<f64>, b: &DVector<f64>, k: usize) -> DVector<f64> {
    let svd = a.clone().svd(true, true);
    let u = svd.u.unwrap();
    let vt = svd.v_t.unwrap();
    let mut x = DVector::zeros(a.ncols());
    for i in 0..k {
        let coef = u.column(i).dot(b) / svd.singular_values[i];   // u_i^T b / sigma_i
        x += coef * vt.row(i).transpose();                       // times v_i
    }
    x
}

fn main() {
    let (s, p_true, b0, b) = signature();
    let svd = s.clone().svd(true, true);
    let u = svd.u.as_ref().unwrap();
    let d = &svd.singular_values;
    for i in 0..5 {
        let beta = u.column(i).dot(&b).abs();
        let beta0 = u.column(i).dot(&b0).abs();
        println!("{:9.4} {:9.4} {:7.4} {:7.4}", d[i], beta, beta / d[i], beta0 / d[i]);
    }
    let x_svd = tsvd(&s, &b, 5);                  // all 5 terms = least squares
    println!("{:e}", (&x_svd - lstsq_qr(&s, &b)).abs().max());
    for k in 1..=5 {
        let xk = tsvd(&s, &b, k);
        let row: Vec<String> = xk.iter().map(|x| format!("{:6.3}", x)).collect();
        println!("{} {}   err {:.3}  res {:.2}", k, row.join(" "),
                 (&xk - &p_true).norm(), (&b - &s * &xk).norm());
    }
}
