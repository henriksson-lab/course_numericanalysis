// C6: CSC arrays with sprs; the inverse of a sparse tridiagonal matrix is dense.
// sprs has no sparse inverse: each column j of A^-1 solves A x = e_j, done here
// with the hand-written Thomas algorithm.
use sprs::{CsMat, TriMat};

/// Thomas algorithm: a = below the diagonal (a[0] unused), b = diagonal,
/// c = above the diagonal (c[n-1] unused), d = right-hand side.
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec()); // copies: inputs stay unchanged
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    x
}

fn main() {
    let mut t = TriMat::new((3, 5));
    for (r, c, v) in [(1, 0, 4.0), (0, 1, 1.0), (2, 1, 2.0), (1, 3, 6.0), (0, 4, 3.0), (2, 4, 1.0)] {
        t.add_triplet(r, c, v);
    }
    let m: CsMat<f64> = t.to_csc();
    for row in m.to_dense().rows() {
        println!("{:?}", row.to_vec());
    }
    println!("{:?} {:?} {:?}", m.data(), m.indices(), m.indptr().raw_storage());
    let lib: Vec<f64> = m.outer_iterator().map(|col| col.data().iter().fold(0.0, |acc, v| acc + v)).collect();
    println!("{:?}", lib);                               // library size of each cell

    let (n, k) = (100, 0.01);
    let mut a = vec![-1.0; n]; a[0] = 0.0;
    let mut c = vec![-1.0; n]; c[n - 1] = 0.0;
    let mut b = vec![2.0 + k; n]; b[0] = 1.0 + k; b[n - 1] = 1.0 + k;
    let nnz_a = 3 * n - 2;
    let (mut nnz_inv, mut min_inv) = (0, f64::INFINITY);
    for j in 0..n {
        let mut e = vec![0.0; n];
        e[j] = 1.0;
        for v in thomas(&a, &b, &c, &e) {
            if v != 0.0 { nnz_inv += 1; }
            min_inv = min_inv.min(v);
        }
    }
    println!("{} {}", nnz_a, nnz_inv);
    println!("{}", min_inv);                             // smallest entry: still positive
    println!("{}", 1e6_f64.powi(2) * 8.0 / 1e12);        // TB for the inverse at 1e6 cells
}
