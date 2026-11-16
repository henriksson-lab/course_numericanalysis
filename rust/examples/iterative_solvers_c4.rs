// rust/examples/iterative_solvers_c4.rs  —  cargo run --release --example iterative_solvers_c4
// C4: CG, Jacobi-preconditioned CG and IC(0)-preconditioned CG on the uniform
// chamber and on the chamber with a cell aggregate (D = 0.1 inside).
include!("iterative_solvers_shared/defs.rs");

fn lap_d(m: usize, dfun: impl Fn(f64, f64) -> f64) -> CsMat<f64> {   // -div(D grad u)
    let (n, h) = (m * m, 1.0 / (m + 1) as f64);
    let mut t = TriMat::new((n, n));
    for j in 0..m {
        for i in 0..m {
            let k = i + m * j;
            let (x, y) = ((i + 1) as f64 * h, (j + 1) as f64 * h);     // node
            let (xw, xe) = ((i as f64 + 0.5) * h, (i as f64 + 1.5) * h); // its faces
            let (ys, yn) = ((j as f64 + 0.5) * h, (j as f64 + 1.5) * h);
            let faces = [(dfun(xw, y), i > 0, k.wrapping_sub(1)),
                         (dfun(xe, y), i < m - 1, k + 1),
                         (dfun(x, ys), j > 0, k.wrapping_sub(m)),
                         (dfun(x, yn), j < m - 1, k + m)];
            for (d, inside, l) in faces {
                t.add_triplet(k, k, d);              // every face adds D to the diagonal
                if inside { t.add_triplet(k, l, -d); }
            }
        }
    }
    t.to_csr()
}
fn aggregate_d(x: f64, y: f64) -> f64 {
    if (x - 0.5).powi(2) + (y - 0.5).powi(2) < 0.3f64.powi(2) { 0.1 } else { 1.0 }
}
// IC(0) for the five-point pattern: diagonal d, west and south entries of L
fn ic0(a: &CsMat<f64>, m: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = a.rows();
    let (mut d, mut lw, mut ls) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let g = |i: usize, j: usize| *a.get(i, j).unwrap_or(&0.0);
    for k in 0..n {
        if k > 0 { lw[k] = g(k, k - 1) / d[k - 1]; }
        if k >= m { ls[k] = g(k, k - m) / d[k - m]; }
        d[k] = (g(k, k) - lw[k].powi(2) - ls[k].powi(2)).sqrt();
    }
    (d, lw, ls)
}
// z = (L L^T)^{-1} r by a forward and a backward substitution
fn ic_solve(l: &(Vec<f64>, Vec<f64>, Vec<f64>), m: usize, r: &[f64]) -> Vec<f64> {
    let (d, lw, ls) = l;
    let n = r.len();
    let mut y = r.to_vec();
    for k in 0..n {
        if k > 0 { y[k] -= lw[k] * y[k - 1]; }
        if k >= m { y[k] -= ls[k] * y[k - m]; }
        y[k] /= d[k];
    }
    for k in (0..n).rev() {
        if k + 1 < n { y[k] -= lw[k + 1] * y[k + 1]; }
        if k + m < n { y[k] -= ls[k + m] * y[k + m]; }
        y[k] /= d[k];
    }
    y
}

fn main() {
    let m = 63;
    let b = chamber_rhs(m, 12.0);
    for a in [lap2(m), lap_d(m, aggregate_d)] {
        let d = diagonal(&a);
        let l = ic0(&a, m);
        let it = |minv: &dyn Fn(&[f64]) -> Vec<f64>| pcg(&a, &b, minv, 1e-6, 10_000).1.unwrap();
        let dmin = d.iter().cloned().fold(f64::INFINITY, f64::min);
        let dmax = d.iter().cloned().fold(0.0, f64::max);
        println!("{} {} {} {} {}", dmin, dmax,
                 it(&|r| r.to_vec()),
                 it(&|r| r.iter().zip(&d).map(|(p, q)| p / q).collect()),
                 it(&|r| ic_solve(&l, m, r)));
    }
}
