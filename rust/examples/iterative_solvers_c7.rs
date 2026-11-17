// rust/examples/iterative_solvers_c7.rs  —  cargo run --release --example iterative_solvers_c7
// C7: hand-written BiCGSTAB (short recurrences, two products with A per step)
// against full GMRES on the perfused chamber with feedback sigma = 0, 100, 200.
include!("iterative_solvers_shared/defs.rs");
use nalgebra::{DMatrix, DVector};

fn adv2(m: usize, pe: f64, sigma: f64) -> CsMat<f64> {   // flow in +x, upwind; minus sigma h^2 I
    let (n, h) = (m * m, 1.0 / (m + 1) as f64);
    let ph = pe * h;
    let mut t = TriMat::new((n, n));
    for j in 0..m {
        for i in 0..m {
            let k = i + m * j;
            t.add_triplet(k, k, 4.0 + ph - sigma * h * h);
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
fn bicgstab(a: &CsMat<f64>, b: &[f64], tol: f64, maxit: usize) -> (Vec<f64>, Option<usize>) {
    let n = b.len();
    let (mut x, mut r, rh) = (vec![0.0; n], b.to_vec(), b.to_vec());
    let nb = norm(b);
    let (mut rho, mut alpha, mut omega) = (1.0, 1.0, 1.0);
    let (mut v, mut p) = (vec![0.0; n], vec![0.0; n]);
    for it in 1..=maxit {
        let rho_new = dot(&rh, &r);
        let beta = (rho_new / rho) * (alpha / omega);
        rho = rho_new;
        for i in 0..n { p[i] = r[i] + beta * (p[i] - omega * v[i]); }
        v = matvec(a, &p);
        alpha = rho / dot(&rh, &v);
        let s: Vec<f64> = (0..n).map(|i| r[i] - alpha * v[i]).collect();
        let t = matvec(a, &s);
        omega = dot(&t, &s) / dot(&t, &t);                 // a GMRES(1) step
        for i in 0..n {
            x[i] += alpha * p[i] + omega * s[i];
            r[i] = s[i] - omega * t[i];
        }
        if norm(&r) < tol * nb { return (x, Some(it)); }
    }
    (x, None)
}
// full GMRES (no restart), as in C5; returns the number of iterations
fn gmres_its(a: &CsMat<f64>, b: &[f64], k: usize, tol: f64) -> usize {
    let beta = norm(b);
    let mut q: Vec<Vec<f64>> = vec![b.iter().map(|v| v / beta).collect()];
    let mut h = DMatrix::<f64>::zeros(k + 1, k);
    for j in 0..k {
        let mut w = matvec(a, &q[j]);
        for i in 0..=j {
            h[(i, j)] = dot(&q[i], &w);
            for (wl, ql) in w.iter_mut().zip(&q[i]) { *wl -= h[(i, j)] * ql; }
        }
        h[(j + 1, j)] = norm(&w);
        q.push(w.iter().map(|v| v / h[(j + 1, j)]).collect());
        let hj = h.view((0, 0), (j + 2, j + 1)).into_owned();
        let mut e1 = DVector::<f64>::zeros(j + 2);
        e1[0] = beta;
        let qr = hj.clone().qr();
        let y = qr.r().solve_upper_triangular(&(qr.q().transpose() * &e1)).unwrap();
        if (&e1 - &hj * &y).norm() < tol * beta { return j + 1; }
    }
    k
}

fn main() {
    let m = 63;
    let b = chamber_rhs(m, 12.0);
    for sigma in [0.0, 100.0, 200.0] {
        let a = adv2(m, 50.0, sigma);
        let (x, it) = bicgstab(&a, &b, 1e-6, 5000);
        let its = match it { Some(i) => format!("{} {}", i, 2 * i), None => "NA NA".to_string() };
        println!("{} {} {:.2e} {}", sigma, its, relres(&a, &b, &x), gmres_its(&a, &b, 300, 1e-6));
    }
}
