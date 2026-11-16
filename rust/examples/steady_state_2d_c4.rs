// C4: assemble the 2D five-point matrix from triplets (k = i + nx * j, x fastest),
// put the edge values into the right-hand side, and check the matrix with a
// manufactured field. sprs has no direct solver: we use a hand-written
// conjugate gradient (CG), which works because the matrix is SPD.
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

fn main() {
    // (a) 3 x 2 grid, h = D = 1, c = 1 on all edges
    let b = edge_rhs(3, 2, 1.0, 1.0, 1.0);
    println!("{:?}", b);
    let (x, _) = cg(&lap2d(3, 2, 1.0, 1.0), &b, 1e-12, 100);
    println!("{:?}", x.iter().map(|v| format!("{:.6}", v)).collect::<Vec<_>>());
    // (b) the chamber
    let (d, lx, ly) = (2000.0_f64, 4000.0_f64, 2000.0_f64);
    let h = 50.0;
    let (nx, ny) = ((lx / h) as usize - 1, (ly / h) as usize - 1);
    let a = lap2d(nx, ny, h, d);
    let sym = (0..nx * ny).all(|k| a.outer_view(k).unwrap().iter()
        .all(|(j, v)| a.get(j, k) == Some(v)));
    println!("{:?} {} {}", a.shape(), a.nnz(), sym);
    // (c) manufactured field, three grids
    let pi = std::f64::consts::PI;
    for h in [100.0_f64, 50.0, 25.0] {
        let (nx, ny) = ((lx / h) as usize - 1, (ly / h) as usize - 1);
        let mut u = vec![0.0; nx * ny];
        for j in 0..ny { for i in 0..nx {
            let (x, y) = (h * (i + 1) as f64, h * (j + 1) as f64);
            u[i + nx * j] = (pi * x / lx).sin() * (2.0 * pi * y / ly).sin();
        } }
        let fac = d * pi * pi * (1.0 / (lx * lx) + 4.0 / (ly * ly));
        let au = matvec(&lap2d(nx, ny, h, d), &u);
        let e = (0..nx * ny).map(|k| (au[k] - fac * u[k]).abs()).fold(0.0, f64::max)
            / u.iter().map(|v| (fac * v).abs()).fold(0.0, f64::max);
        println!("h = {:3}: relative error {:.3e}", h, e);
    }
}
