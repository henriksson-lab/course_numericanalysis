// C5: ridge regularisation two ways, and the numbers behind the L-curve.
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

fn main() {
    let (s, _p_true, _b0, b) = signature();
    let x1 = ridge(&s, &b, 10.0);
    // stacked least squares: [S; lambda I] x ~ [b; 0]
    let mut a_aug = DMatrix::zeros(45, 5);
    a_aug.view_mut((0, 0), (40, 5)).copy_from(&s);
    a_aug.view_mut((40, 0), (5, 5)).copy_from(&(10.0 * DMatrix::<f64>::identity(5, 5)));
    let mut b_aug = DVector::zeros(45);
    b_aug.rows_mut(0, 40).copy_from(&b);
    let x2 = lstsq_qr(&a_aug, &b_aug);
    println!("{}", fmt(&x1, 4));
    println!("{:e}", (&x1 - &x2).abs().max());
    for k in 0..9 {
        let lambda = 10f64.powf(-1.0 + 0.5 * k as f64);
        let x = ridge(&s, &b, lambda);
        let row: Vec<String> = x.iter().map(|v| format!("{:6.3}", v)).collect();
        println!("{:8.2} {}  res {:6.2}  norm {:.3}  sum {:.3}", lambda, row.join(" "),
                 (&b - &s * &x).norm(), x.norm(), x.sum());
    }
}
