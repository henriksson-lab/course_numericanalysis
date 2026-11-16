// Background: Cholesky with nalgebra, and a small sparse matrix with sprs.
use nalgebra::{dmatrix, dvector};
use sprs::TriMat;

fn main() {
    let a = dmatrix![4.0, 2.0; 2.0, 5.0];
    let chol = a.cholesky().expect("not positive definite"); // None if not SPD
    println!("{}", chol.l());          // lower triangular: A = L L^T
    let x = chol.solve(&dvector![2.0, 9.0]);
    println!("{:?}", x.as_slice());

    // sprs builds CSC from triplets (row, column, value), 0-based
    let mut t = TriMat::new((4, 4));
    for (r, c, v) in [(0, 0, 3.0), (2, 0, 1.0), (1, 1, 5.0), (0, 3, 2.0), (3, 3, 7.0)] {
        t.add_triplet(r, c, v);
    }
    let m: sprs::CsMat<f64> = t.to_csc();
    println!("{}", m.nnz());
}
