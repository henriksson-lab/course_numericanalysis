// C3: Lauchli's example. The normal equations lose delta^2; QR does not.
use nalgebra::{DMatrix, DVector};

fn main() {
    for d in [1e-6, 1e-8] {
        let a = DMatrix::from_row_slice(3, 2, &[1.0, 1.0, d, 0.0, 0.0, d]);
        let b = DVector::from_vec(vec![2.0, d, d]);          // exact solution x = (1, 1)
        match (a.transpose() * &a).lu().solve(&(a.transpose() * &b)) {
            Some(x) => println!("{:.15} {:.15}", x[0], x[1]),
            None => println!("normal equations: singular"),
        }
        let qr = a.clone().qr();
        let mut qtb = b.clone();
        qr.q_tr_mul(&mut qtb);
        let x = qr.r().solve_upper_triangular(&qtb.rows(0, 2).into_owned()).unwrap();
        println!("{:.15} {:.15}", x[0], x[1]);
    }
}
