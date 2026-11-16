// C5: the ordering bug. The matrix numbers the grid with x fastest (k = i + nx * j),
// but the vectors are filled with y fastest (k = j + ny * i), as numpy's default
// ravel would. lap2d, edge_rhs, matvec, dot and cg are as in steady_state_2d_c4.rs.
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

fn edge_rhs(nx: usize, ny: usize, h: f64, d: f64, c0: f64, wrong: bool) -> Vec<f64> {
    let s = d / (h * h);
    let mut b = vec![0.0; nx * ny];
    for j in 0..ny {
        for i in 0..nx {
            let mut m = 0.0;                       // number of edge neighbours
            if i == 0 { m += 1.0; }
            if i == nx - 1 { m += 1.0; }
            if j == 0 { m += 1.0; }
            if j == ny - 1 { m += 1.0; }
            let k = if wrong { j + ny * i } else { i + nx * j };
            b[k] = s * c0 * m;
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

fn check(lx: f64, ly: f64, h: f64, my: f64, wrong: bool) {
    let (d, c0) = (2000.0, 200.0);
    let pi = std::f64::consts::PI;
    let (nx, ny) = ((lx / h) as usize - 1, (ly / h) as usize - 1);
    let a = lap2d(nx, ny, h, d);                       // x fastest
    let mut u = vec![0.0; nx * ny];
    for j in 0..ny { for i in 0..nx {
        let (x, y) = (h * (i + 1) as f64, h * (j + 1) as f64);
        let k = if wrong { j + ny * i } else { i + nx * j };
        u[k] = (pi * x / lx).sin() * (my * pi * y / ly).sin();
    } }
    let fac = d * pi * pi * (1.0 / (lx * lx) + my * my / (ly * ly));
    let au = matvec(&a, &u);
    let err = (0..nx * ny).map(|k| (au[k] - fac * u[k]).abs()).fold(0.0, f64::max)
        / u.iter().map(|v| (fac * v).abs()).fold(0.0, f64::max);
    let b: Vec<f64> = edge_rhs(nx, ny, h, d, c0, wrong).iter().map(|v| v - 0.5).collect();
    let (c, _) = cg(&a, &b, 1e-12, 10000);
    let cmax = c.iter().cloned().fold(f64::MIN, f64::max);
    println!("{:4} x {:4}, my = {}, wrong = {:>5}: check {:.3e}, max c {:.1}",
             lx, ly, my, wrong, err, cmax);
}

fn main() {
    check(4000.0, 2000.0, 50.0, 2.0, false);
    check(4000.0, 2000.0, 50.0, 2.0, true);
    check(2000.0, 2000.0, 50.0, 2.0, true);            // square chamber: the bug is invisible
}
