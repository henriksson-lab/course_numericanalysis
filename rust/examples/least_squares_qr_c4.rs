// C4: classical and modified Gram-Schmidt, and loss of orthogonality.
use nalgebra::{DMatrix, DVector};

// scientific notation as R and Python print it: 2.2e+01 rather than Rust's 2.2e1
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn gs(a: &DMatrix<f64>, modified: bool) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = a.ncols();
    let mut q = a.clone();
    let mut r = DMatrix::zeros(n, n);
    for j in 0..n {
        let mut v: DVector<f64> = a.column(j).into_owned();
        for i in 0..j {
            r[(i, j)] = if modified { q.column(i).dot(&v) } else { q.column(i).dot(&a.column(j)) };
            v -= r[(i, j)] * q.column(i);
        }
        r[(j, j)] = v.norm();
        q.set_column(j, &(v / r[(j, j)]));
    }
    (q, r)
}

fn orth(q: &DMatrix<f64>) -> f64 {
    // 2-norm of Q^T Q - I = its largest singular value
    let e = q.transpose() * q - DMatrix::identity(q.ncols(), q.ncols());
    e.singular_values().max()
}

fn kappa(a: &DMatrix<f64>) -> f64 {
    let sv = a.clone().singular_values();
    sv.max() / sv.min()
}

fn main() {
    // B3: columns 1 and (1, 3, 5)
    let a3 = DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 1.0, 3.0, 1.0, 5.0]);
    let (q, r) = gs(&a3, false);
    println!("{:.6} {:.6}\n{:.6} {:.6}", r[(0, 0)], r[(0, 1)], r[(1, 0)], r[(1, 1)]);
    let qtb = q.transpose() * DVector::from_vec(vec![35.0, 28.4, 21.6]);
    let x = r.solve_upper_triangular(&qtb).unwrap();
    println!("{:.5} {:.5}", x[0], x[1]);

    let t: Vec<f64> = (0..50).map(|i| i as f64 / 49.0).collect();
    for d in [2, 4, 6, 8, 10, 12] {
        let a = DMatrix::from_fn(50, d + 1, |i, j| t[i].powi(j as i32));
        let q_h = a.clone().qr().q();                          // Householder
        println!("{:2} {:>8} {:>8} {:>8} {:>8}", d, sci(kappa(&a), 1),
                 sci(orth(&gs(&a, false).0), 1), sci(orth(&gs(&a, true).0), 1), sci(orth(&q_h), 1));
    }
}
