// C6: residual and error of LU solve vs multiplying by the inverse.
use nalgebra::{dmatrix, dvector};

fn main() {
    let f = dvector![0.3, 0.2, 0.4, 0.1];
    for p in 0..=6 {
        let d = 10f64.powi(-p);
        let a4 = dmatrix![  1.0,  1.0,           1.0,   1.0;
                          -26.0, -26.0 + d,    -10.0, -18.0;   // d13C: sources 1, 2 differ by d
                            2.0,   2.0 + d,      6.0,  14.0;   // d15N
                            4.0,   4.0 + 2.0 * d, 3.0, 19.0];  // d34S
        let b = &a4 * &f;
        let lu = a4.clone().lu();
        let x_lu = lu.solve(&b).unwrap();
        let x_inv = lu.try_inverse().unwrap() * &b;   // inverse via the same LU
        println!("d = {:.0e}  residual LU {:.1e}  inv {:.1e}   error LU {:.1e}  inv {:.1e}", d,
                 (&b - &a4 * &x_lu).abs().max(), (&b - &a4 * &x_inv).abs().max(),
                 (&x_lu - &f).abs().max(), (&x_inv - &f).abs().max());
    }
}
