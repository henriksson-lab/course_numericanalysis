// rust/examples/iterative_solvers_c6.rs  —  cargo run --release --example iterative_solvers_c6
// C6: restarted GMRES(k) on the perfused chamber and on the signal chamber
// (positive feedback, sigma = 100), where short restarts stall. Hand-written
// GMRES as in C5, plus a restart loop; the smooth bump v tests v'Av > 0.
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
fn residual(a: &CsMat<f64>, b: &[f64], x: &[f64]) -> Vec<f64> {
    let ax = matvec(a, x);
    b.iter().zip(&ax).map(|(p, q)| p - q).collect()
}
// one GMRES cycle of at most k steps from x = 0; stops when ||b - A x|| < tol ||b||
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
// restarted GMRES(k): at most maxit steps in total; returns iterations and relative residual
fn gmres_k(a: &CsMat<f64>, b: &[f64], k: usize, tol: f64, maxit: usize) -> (usize, f64) {
    let (mut x, nb, mut it) = (vec![0.0; b.len()], norm(b), 0);
    loop {
        let r = residual(a, b, &x);
        let nr = norm(&r);
        if nr < tol * nb || it >= maxit { return (it, nr / nb); }
        let (dx, j) = gmres(a, &r, k.min(maxit - it), tol * nb / nr);   // one cycle
        for (xi, d) in x.iter_mut().zip(&dx) { *xi += d; }
        it += j;
    }
}
fn sig3(v: f64) -> f64 {                                // round to 3 significant digits
    let e = 10f64.powi(2 - v.abs().log10().floor() as i32);
    (v * e).round() / e
}

fn main() {
    let m = 63;
    let h = 1.0 / (m + 1) as f64;
    let b = chamber_rhs(m, 12.0);
    let a = adv2(m, 50.0, 0.0);                         // perfused chamber
    let s = adv2(m, 50.0, 100.0);                       // signal with positive feedback
    for k in [10, 20, 50, 300] {
        let (it1, _) = gmres_k(&a, &b, k, 1e-6, 3000);
        let (it2, r2) = gmres_k(&s, &b, k, 1e-6, 3000);
        if r2 > 1e-3 { println!("{} {} {} {}", k, it1, it2, sig3(r2)); }
        else { println!("{} {} {} {:e}", k, it1, it2, sig3(r2)); }
    }
    let sv: Vec<f64> = (1..=m).map(|i| (std::f64::consts::PI * i as f64 * h).sin()).collect();
    let v: Vec<f64> = (0..m * m).map(|k| sv[k % m] * sv[k / m]).collect();   // smooth bump
    let vv = dot(&v, &v);
    println!("{} {}", dot(&v, &matvec(&a, &v)) / vv, dot(&v, &matvec(&s, &v)) / vv);
}
