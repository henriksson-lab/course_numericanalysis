// The SVD of a 2x2 matrix with nalgebra (R: svd; numpy: np.linalg.svd).
use nalgebra::{DMatrix, DVector};

fn main() {
    let a = DMatrix::from_row_slice(2, 2, &[3.0, 0.0,
                                            4.0, 5.0]);
    let svd = a.clone().svd(true, true);     // sorted, largest first
    let d: DVector<f64> = svd.singular_values.clone();
    let u = svd.u.unwrap();
    let vt = svd.v_t.unwrap();               // like numpy, nalgebra returns V^T
    println!("{} {}", d[0], d[1]);           // singular values
    println!("{}", d[0] / d[1]);             // condition number
    let rebuilt = &u * DMatrix::from_diagonal(&d) * &vt;
    println!("{:e}", (a - rebuilt).abs().max()); // A = U Sigma V^T
}
