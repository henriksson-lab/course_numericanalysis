// rust/examples/iterative_solvers_c5.rs  —  cargo run --release --example iterative_solvers_c5
// C5: the perfused chamber (advection-diffusion, nonsymmetric). CG fails;
// hand-written GMRES (Arnoldi on sprs, small least squares by nalgebra QR) works.
include!("iterative_solvers_shared/defs.rs");
use nalgebra::{DMatrix, DVector};

fn adv2(m: usize, pe: f64) -> CsMat<f64> {          // flow in +x, upwind differences
    let (n, ph) = (m * m, pe / (m + 1) as f64);
    let mut t = TriMat::new((n, n));
    for j in 0..m {
        for i in 0..m {
            let k = i + m * j;
            t.add_triplet(k, k, 4.0 + ph);
            if i > 0 { t.add_triplet(k, k - 1, -1.0 - ph); }
            if i < m - 1 { t.add_triplet(k, k + 1, -1.0); }
            if j > 0 { t.add_triplet(k, k - m, -1.0); }
            if j < m - 1 { t.add_triplet(k, k + m, -1.0); }
        }
    }
    t.to_csr()
}
fn relres(a: &CsMat<f64>, b: &[f64], x: &[f64]) -> f64 {
    let ax = matvec(a, x);
    norm(&b.iter().zip(&ax).map(|(p, q)| p - q).collect::<Vec<_>>()) / norm(b)
}
fn cg_fixed(a: &CsMat<f64>, b: &[f64], steps: usize) -> Vec<f64> {   // plain CG, no stopping
    let n = b.len();
    let (mut x, mut r, mut p) = (vec![0.0; n], b.to_vec(), b.to_vec());
    let mut rr = dot(&r, &r);
    for _ in 0..steps {
        let ap = matvec(a, &p);
        let alpha = rr / dot(&p, &ap);
        for i in 0..n { x[i] += alpha * p[i]; r[i] -= alpha * ap[i]; }
        let rr_new = dot(&r, &r);
        for i in 0..n { p[i] = r[i] + (rr_new / rr) * p[i]; }
        rr = rr_new;
    }
    x
}
fn gmres(a: &CsMat<f64>, b: &[f64], k: usize, tol: f64) -> (Vec<f64>, usize) {
    let beta = norm(b);
    let mut q: Vec<Vec<f64>> = vec![b.iter().map(|v| v / beta).collect()];
    let mut h = DMatrix::<f64>::zeros(k + 1, k);
    let mut y = DVector::<f64>::zeros(0);
    let mut jlast = 0;
    for j in 0..k {
        let mut w = matvec(a, &q[j]);                   // one product with A
        for i in 0..=j {
            h[(i, j)] = dot(&q[i], &w);
            for (wl, ql) in w.iter_mut().zip(&q[i]) { *wl -= h[(i, j)] * ql; }
        }
        h[(j + 1, j)] = norm(&w);
        q.push(w.iter().map(|v| v / h[(j + 1, j)]).collect());
        // small least-squares problem: minimise || beta e1 - H y ||
        let hj = h.view((0, 0), (j + 2, j + 1)).into_owned();
        let mut e1 = DVector::<f64>::zeros(j + 2);
        e1[0] = beta;
        let qr = hj.clone().qr();
        y = qr.r().solve_upper_triangular(&(qr.q().transpose() * &e1)).unwrap();
        jlast = j + 1;
        if (&e1 - &hj * &y).norm() < tol * beta { break; }
    }
    let mut x = vec![0.0; b.len()];
    for (l, yl) in y.iter().enumerate() {
        for (xi, qi) in x.iter_mut().zip(&q[l]) { *xi += yl * qi; }
    }
    (x, jlast)
}

fn main() {
    let m = 63;
    let b = chamber_rhs(m, 12.0);
    let a = adv2(m, 50.0);
    let at: CsMat<f64> = a.transpose_view().to_csr();
    println!("{}", a == at);                            // symmetric?
    let res: Vec<f64> = [50, 100, 300].iter().map(|&s| relres(&a, &b, &cg_fixed(&a, &b, s))).collect();
    println!("{:?}", res);
    let (x, it) = gmres(&a, &b, 300, 1e-6);
    println!("{} {:e}", it, relres(&a, &b, &x));
    println!("{}", gmres(&lap2(m), &b, 300, 1e-6).1);  // GMRES on the symmetric chamber
    println!("{}", 1.0 - x[(m * m - 1) / 2]);
}
