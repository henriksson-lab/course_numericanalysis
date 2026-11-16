// C6: oxygen under a lid with Michaelis-Menten uptake, solved by Newton.
// Each Newton step solves J delta = -F with J = A + diag(r'(c)), which is SPD,
// so the hand-written CG works (sprs has no direct solver).
// lap2d, edge_rhs, matvec, dot and cg are as in steady_state_2d_c4.rs.
use sprs::{CsMat, TriMat};

fn lap2d(nx: usize, ny: usize, h: f64, d: f64) -> CsMat<f64> {
    let s = d / (h * h);
    let idx = |i: usize, j: usize| i + nx * j;
    let mut t = TriMat::new((nx * ny, nx * ny));
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            t.add_triplet(k, k, 4.0 * s);
            if i > 0 { t.add_triplet(k, idx(i - 1, j), -s); }
            if i + 1 < nx { t.add_triplet(k, idx(i + 1, j), -s); }
            if j > 0 { t.add_triplet(k, idx(i, j - 1), -s); }
            if j + 1 < ny { t.add_triplet(k, idx(i, j + 1), -s); }
        }
    }
    t.to_csr()
}

fn edge_rhs(nx: usize, ny: usize, h: f64, d: f64, c0: f64) -> Vec<f64> {
    let s = d / (h * h);
    let mut b = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let mut m = 0.0;                       // number of edge neighbours
            if i == 0 { m += 1.0; }
            if i == nx - 1 { m += 1.0; }
            if j == 0 { m += 1.0; }
            if j == ny - 1 { m += 1.0; }
            b[i + nx * j] = s * c0 * m;
        }
    }
    b
}

fn matvec(a: &CsMat<f64>, x: &[f64]) -> Vec<f64> {
    a.outer_iterator().map(|row| row.iter().map(|(j, v)| v * x[j]).sum()).collect()
}
fn dot(u: &[f64], v: &[f64]) -> f64 { u.iter().zip(v).map(|(p, q)| p * q).sum() }

// conjugate gradient for SPD a; stops when ||r|| < tol * ||b||
fn cg(a: &CsMat<f64>, b: &[f64], tol: f64, maxit: usize) -> (Vec<f64>, usize) {
    let n = b.len();
    let (mut x, mut r, mut p) = (vec![0.0; n], b.to_vec(), b.to_vec());
    let mut rr = dot(&r, &r);
    let stop = tol * rr.sqrt();
    for it in 1..=maxit {
        let ap = matvec(a, &p);
        let alpha = rr / dot(&p, &ap);
        for i in 0..n { x[i] += alpha * p[i]; r[i] -= alpha * ap[i]; }
        let rr_new = dot(&r, &r);
        if rr_new.sqrt() < stop { return (x, it); }
        for i in 0..n { p[i] = r[i] + (rr_new / rr) * p[i]; }
        rr = rr_new;
    }
    panic!("CG did not converge in {} iterations", maxit);
}

// returns (c with k = i + nx * j, Newton steps)
fn oxygen_mm(a: &CsMat<f64>, b: &[f64], r0: f64, km: f64, c0: f64, clip: bool) -> (Vec<f64>, usize) {
    let n = b.len();
    let mut c = vec![c0; n];
    let mut it = 0;
    loop {
        let ac = matvec(a, &c);
        let f: Vec<f64> = (0..n).map(|k| ac[k] + r0 * c[k] / (km + c[k]) - b[k]).collect();
        if dot(&f, &f).sqrt() < 1e-8 || it == 50 { break; }
        // J = A + diag(r'(c)): add the diagonal to a copy of A
        let mut t = TriMat::new((n, n));
        for (v, (i, j)) in a.iter() { t.add_triplet(i, j, *v); }
        for k in 0..n { t.add_triplet(k, k, r0 * km / ((km + c[k]) * (km + c[k]))); }
        let jac: CsMat<f64> = t.to_csr();
        let minus_f: Vec<f64> = f.iter().map(|v| -v).collect();
        let (delta, _) = cg(&jac, &minus_f, 1e-12, 20000);
        for k in 0..n {
            c[k] += delta[k];
            if clip { c[k] = c[k].max(0.0); }
        }
        it += 1;
    }
    (c, it)
}

fn main() {
    let (d, c0, km, lx, ly, h) = (2000.0_f64, 200.0, 1.0, 4000.0_f64, 2000.0_f64, 50.0);
    let (nx, ny) = ((lx / h) as usize - 1, (ly / h) as usize - 1);
    let a = lap2d(nx, ny, h, d);
    let b = edge_rhs(nx, ny, h, d, c0);
    let min = |c: &[f64]| c.iter().cloned().fold(f64::MAX, f64::min);
    let (c, it) = oxygen_mm(&a, &b, 1.0, km, c0, false);
    println!("plain Newton, 2e5 cells/cm2: {} steps, min c = {:.2}", it, min(&c));
    for rho in [1.0, 2.0, 3.0] {
        let (c, it) = oxygen_mm(&a, &b, 0.5 * rho, km, c0, true);
        let hyp = c.iter().filter(|&&v| v < 10.0).count() as f64 / c.len() as f64;
        println!("{} e5 cells/cm2: {:2} steps, min c = {:7.3} uM, hypoxic {:4.1} %",
                 rho, it, min(&c), 100.0 * hyp);
    }
}
