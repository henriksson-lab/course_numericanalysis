// Shared by the iterative_solvers_* examples: include!("iterative_solvers_shared/defs.rs")
// The oxygen-in-a-chamber matrix (2D five-point Laplacian, natural ordering
// k = i + m*j), its right-hand side, and small vector helpers.
use sprs::{CsMat, TriMat};

#[allow(dead_code)]
fn lap2(m: usize) -> CsMat<f64> {                 // 4 on the diagonal, -1 for neighbours
    let n = m * m;
    let mut t = TriMat::new((n, n));
    for j in 0..m {
        for i in 0..m {
            let k = i + m * j;
            t.add_triplet(k, k, 4.0);
            if i > 0 { t.add_triplet(k, k - 1, -1.0); }
            if i < m - 1 { t.add_triplet(k, k + 1, -1.0); }
            if j > 0 { t.add_triplet(k, k - m, -1.0); }
            if j < m - 1 { t.add_triplet(k, k + m, -1.0); }
        }
    }
    t.to_csr()
}
#[allow(dead_code)]
fn chamber_rhs(m: usize, phi: f64) -> Vec<f64> {
    vec![phi / ((m + 1) * (m + 1)) as f64; m * m]
}
#[allow(dead_code)]
fn matvec(a: &CsMat<f64>, x: &[f64]) -> Vec<f64> {   // CSR: one row at a time
    a.outer_iterator().map(|row| row.iter().map(|(j, v)| v * x[j]).sum()).collect()
}
#[allow(dead_code)]
fn dot(u: &[f64], v: &[f64]) -> f64 { u.iter().zip(v).map(|(p, q)| p * q).sum() }
#[allow(dead_code)]
fn norm(v: &[f64]) -> f64 { dot(v, v).sqrt() }
#[allow(dead_code)]
fn diagonal(a: &CsMat<f64>) -> Vec<f64> {
    (0..a.rows()).map(|i| *a.get(i, i).unwrap_or(&0.0)).collect()
}
// Preconditioned CG; minv(r) applies M^{-1}. Stops when ||r|| < tol ||b||.
#[allow(dead_code)]
fn pcg(a: &CsMat<f64>, b: &[f64], minv: impl Fn(&[f64]) -> Vec<f64>,
       tol: f64, maxit: usize) -> (Vec<f64>, Option<usize>) {
    let n = b.len();
    let (mut x, mut r) = (vec![0.0; n], b.to_vec());
    let mut z = minv(&r);
    let mut p = z.clone();
    let (mut rz, nb) = (dot(&r, &z), norm(b));
    for it in 1..=maxit {
        let ap = matvec(a, &p);
        let alpha = rz / dot(&p, &ap);
        for i in 0..n { x[i] += alpha * p[i]; r[i] -= alpha * ap[i]; }
        if norm(&r) < tol * nb { return (x, Some(it)); }
        z = minv(&r);
        let rz_new = dot(&r, &z);
        for i in 0..n { p[i] = z[i] + (rz_new / rz) * p[i]; }
        rz = rz_new;
    }
    (x, None)
}
