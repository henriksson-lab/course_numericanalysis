// C2: least squares fails on the near rank-deficient signature; the SVD shows why.
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

fn main() {
    let (s, p_true, _b0, b) = signature();
    let x = lstsq_qr(&s, &b);
    println!("{}", fmt(&x, 3));
    println!("{}", (&b - &s * &x).norm());
    println!("{}", (&b - &s * &p_true).norm());
    let svd = s.clone().svd(true, true);
    let d = &svd.singular_values;
    println!("{}", fmt(d, 6));
    println!("{}", d[0] / d[4]);
    let vt = svd.v_t.as_ref().unwrap();
    let mut v5: DVector<f64> = vt.row(4).transpose();
    if v5[0] < 0.0 { v5 = -v5; }                 // sign: CD4 entry positive
    println!("{}", fmt(&v5, 4));
    println!("{}", (&s * (0.1 * &v5)).norm());   // a step of 0.1 along v5
    let e_b = DVector::from_vec(vec![0.0, 0.0, 0.1, 0.0, 0.0]);
    println!("{}", (&s * e_b).norm());           // 0.1 more B cells
}
