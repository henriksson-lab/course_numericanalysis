// C7: fill-in for an arrow matrix (hub gene). sprs has no sparse Cholesky
// (the separate sprs-ldl crate has one), so we run nalgebra's dense Cholesky and
// count the non-zeros of L: the fill-in pattern is the same.
use nalgebra::DMatrix;

fn arrow(n: usize, hub: usize) -> DMatrix<f64> {
    DMatrix::from_fn(n, n, |i, j| {
        if i == hub && j == hub { n as f64 }
        else if i == hub || j == hub { 1.0 }
        else if i == j { 2.0 } else { 0.0 }
    })
}

fn nnz(m: &DMatrix<f64>) -> usize {
    m.iter().filter(|v| **v != 0.0).count()
}

fn main() {
    let n = 1000;
    println!("{}", nnz(&arrow(n, 0)));
    let l_first = arrow(n, 0).cholesky().unwrap().l();    // hub first
    let l_last = arrow(n, n - 1).cholesky().unwrap().l(); // hub last
    println!("{} {}", nnz(&l_first), nnz(&l_last));
}
