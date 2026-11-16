// C1: Cholesky of a 3 x 3 covariance matrix, triangular solves, log det,
// and a "correlation matrix" that is not positive definite.
use nalgebra::{dmatrix, dvector};

fn main() {
    let sigma = dmatrix![4.0_f64, 2.0, -2.0; 2.0, 10.0, 2.0; -2.0, 2.0, 6.0];
    let l = sigma.clone().cholesky().expect("not SPD").l();
    println!("{}", l);                                   // compare with B1(a)
    let b = dvector![4.0, 14.0, 6.0];
    let y = l.solve_lower_triangular(&b).unwrap();      // B1(b)
    println!("{:?}", y.as_slice());
    let x = l.transpose().solve_upper_triangular(&y).unwrap();
    println!("{:?}", x.as_slice());
    let logdet: f64 = 2.0 * l.diagonal().iter().map(|v| v.ln()).sum::<f64>();
    println!("{}", logdet);                              // log det, B1(c)
    println!("{}", sigma.determinant().ln());

    let c = dmatrix![1.0_f64, 0.9, 0.9; 0.9, 1.0, -0.9; 0.9, -0.9, 1.0];
    match c.clone().cholesky() {                         // B1(d)
        Some(_) => println!("SPD"),
        None => println!("cholesky: None (not positive definite)"),
    }
    let mut ev: Vec<f64> = c.symmetric_eigen().eigenvalues.iter().cloned().collect();
    ev.sort_by(|p, q| q.partial_cmp(p).unwrap());
    println!("{:?}", ev);
}
