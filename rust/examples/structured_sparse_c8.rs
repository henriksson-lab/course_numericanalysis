// C8: Jacobi and conjugate gradient (hand-written) on the morphogen system,
// with the matrix stored in CSR by sprs and a hand-written matrix-vector product.
use sprs::{CsMat, TriMat};

fn matvec(a: &CsMat<f64>, x: &[f64]) -> Vec<f64> {   // CSR: one row at a time
    a.outer_iterator()
        .map(|row| row.iter().map(|(j, v)| v * x[j]).sum())
        .collect()
}
fn norm(v: &[f64]) -> f64 { v.iter().map(|x| x * x).sum::<f64>().sqrt() }
fn dot(u: &[f64], v: &[f64]) -> f64 { u.iter().zip(v).map(|(p, q)| p * q).sum() }

fn jacobi(a: &CsMat<f64>, b: &[f64], tol: f64, maxit: usize) -> (Vec<f64>, Option<usize>) {
    let n = b.len();
    let d: Vec<f64> = (0..n).map(|i| *a.get(i, i).unwrap()).collect();
    let mut x = vec![0.0; n];
    for it in 1..=maxit {
        let ax = matvec(a, &x);
        for i in 0..n { x[i] += (b[i] - ax[i]) / d[i]; }
        let ax = matvec(a, &x);
        let r: Vec<f64> = (0..n).map(|i| b[i] - ax[i]).collect();
        if norm(&r) < tol { return (x, Some(it)); }
    }
    (x, None)                                           // not converged
}

fn cg(a: &CsMat<f64>, b: &[f64], tol: f64, maxit: usize) -> (Vec<f64>, Option<usize>) {
    let n = b.len();
    let (mut x, mut r, mut p) = (vec![0.0; n], b.to_vec(), b.to_vec());
    let mut rr = dot(&r, &r);
    for it in 1..=maxit {
        let ap = matvec(a, &p);
        let alpha = rr / dot(&p, &ap);
        for i in 0..n { x[i] += alpha * p[i]; r[i] -= alpha * ap[i]; }
        let rr_new = dot(&r, &r);
        if rr_new.sqrt() < tol { return (x, Some(it)); }
        for i in 0..n { p[i] = r[i] + (rr_new / rr) * p[i]; }
        rr = rr_new;
    }
    (x, None)
}

fn main() {
    let (n, k) = (1000, 0.01);
    let mut t = TriMat::new((n, n));
    for i in 0..n {
        let end = i == 0 || i == n - 1;
        t.add_triplet(i, i, if end { 1.0 + k } else { 2.0 + k });
        if i > 0 { t.add_triplet(i, i - 1, -1.0); }
        if i < n - 1 { t.add_triplet(i, i + 1, -1.0); }
    }
    let a: CsMat<f64> = t.to_csr();
    let mut s = vec![0.0; n];
    s[0] = 1.0;
    let (xj, itj) = jacobi(&a, &s, 1e-10, 20000);
    let (xc, itc) = cg(&a, &s, 1e-10, 20000);
    println!("{:?} {:?}", itj, itc);
    let diff = xj.iter().zip(&xc).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
    println!("{:e}", diff);                              // Jacobi vs CG

    let mut tb = TriMat::new((2, 2));
    for (i, j, v) in [(0, 0, 1.0), (0, 1, 2.0), (1, 0, 2.0), (1, 1, 1.0)] { tb.add_triplet(i, j, v); }
    let bm: CsMat<f64> = tb.to_csr();
    let (xb, itb) = jacobi(&bm, &[3.0, 3.0], 1e-10, 50);
    println!("{:?}", itb);                               // does not converge
    println!("{:e}", xb.iter().fold(0.0_f64, |m, v| m.max(v.abs())));
}
