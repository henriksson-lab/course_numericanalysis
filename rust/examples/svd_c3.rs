// C3: exact rank deficiency (CD8 identical to CD4): numerical rank, pseudoinverse.
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

// pseudoinverse by the SVD: drop sigma_i <= rtol * sigma_1
fn pinv(a: &DMatrix<f64>, rtol: f64) -> DMatrix<f64> {
    let svd = a.clone().svd(true, true);
    let d = &svd.singular_values;
    let tol = rtol * d[0];
    let dplus = DMatrix::from_diagonal(&d.map(|s| if s > tol { 1.0 / s } else { 0.0 }));
    svd.v_t.as_ref().unwrap().transpose() * dplus * svd.u.as_ref().unwrap().transpose()
}

fn main() {
    let (s, _p_true, _b0, b) = signature();
    let mut s0 = s.clone();
    s0.set_column(1, &s.column(0));                 // CD8 identical to CD4
    let d0 = s0.clone().singular_values();
    println!("{:?}", d0.as_slice());
    let rtol = 40.0 * f64::EPSILON;
    println!("{}", d0.iter().filter(|&&x| x > rtol * d0[0]).count());   // numerical rank

    let x0 = pinv(&s0, rtol) * &b;                  // minimum-norm solution
    let s0_drop = s0.clone().remove_column(1);      // drop CD8, as R's lm does
    let xd4 = lstsq_qr(&s0_drop, &b);
    let xd = DVector::from_vec(vec![xd4[0], 0.0, xd4[1], xd4[2], xd4[3]]);
    println!("{}", fmt(&x0, 4));
    println!("{}", fmt(&xd, 4));
    println!("{} {}", (&b - &s0 * &x0).norm(), (&b - &s0 * &xd).norm());  // same fit
    println!("{} {}", x0.norm(), xd.norm());                              // different norms
    // nalgebra's own: pseudo_inverse takes an ABSOLUTE tolerance
    let x_na = s0.clone().pseudo_inverse(rtol * d0[0]).unwrap() * &b;
    println!("{}", fmt(&x_na, 4));
}
