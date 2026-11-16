// Q5: which triangle? R's chol() returns the upper factor R = L^T.
use nalgebra::dmatrix;

fn main() {
    let sigma = dmatrix![4.0, 2.0, -2.0; 2.0, 10.0, 2.0; -2.0, 2.0, 6.0];
    let l = sigma.cholesky().unwrap().l();   // nalgebra: lower factor
    let r = l.transpose();                   // what R's chol() returns
    println!("{}", &r * r.transpose());      // covariance of R z: wrong
    println!("{}", &l * l.transpose());      // covariance of L z: Sigma
}
