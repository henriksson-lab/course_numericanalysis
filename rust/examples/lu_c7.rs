// C7: a zero pivot and a tiny pivot break LU without pivoting.
use nalgebra::{dmatrix, dvector, DMatrix, DVector};

fn lu_nopivot(a: &DMatrix<f64>) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = a.nrows();
    let mut l = DMatrix::<f64>::identity(n, n);
    let mut u = a.clone();
    for k in 0..n - 1 {
        for i in k + 1..n {
            l[(i, k)] = u[(i, k)] / u[(k, k)];
            for j in k..n {
                u[(i, j)] -= l[(i, k)] * u[(k, j)];
            }
        }
    }
    (l, u)
}

fn forward_sub(l: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let n = b.len();
    let mut y = DVector::<f64>::zeros(n);
    for i in 0..n {
        let s: f64 = (0..i).map(|j| l[(i, j)] * y[j]).sum();
        y[i] = (b[i] - s) / l[(i, i)];
    }
    y
}

fn back_sub(u: &DMatrix<f64>, y: &DVector<f64>) -> DVector<f64> {
    let n = y.len();
    let mut x = DVector::<f64>::zeros(n);
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|j| u[(i, j)] * x[j]).sum();
        x[i] = (y[i] - s) / u[(i, i)];
    }
    x
}

fn show(m: &DMatrix<f64>) {
    for r in m.row_iter() {
        println!("{:?}", r.iter().collect::<Vec<_>>());
    }
}

fn main() {
    let a0 = dmatrix![1.0, 1.0, 1.0; -26.0, -26.0, -10.0; 2.0, 10.0, 6.0];
    let b0 = dvector![1.0, -20.0, 6.5];
    let (_l0, u0) = lu_nopivot(&a0);
    println!("{}", u0.iter().all(|v| v.is_finite()));      // the check
    let mut a0s = a0.clone();
    a0s.swap_rows(1, 2);                                   // swap rows 2 and 3
    let mut b0s = b0.clone();
    b0s.swap_rows(1, 2);
    let (l0s, u0s) = lu_nopivot(&a0s);
    show(&l0s);
    show(&u0s);
    println!("{:?}", back_sub(&u0s, &forward_sub(&l0s, &b0s)).as_slice());

    let a_s = dmatrix![1e-20, 1.0; 1.0, 1.0];
    let b_s = dvector![1.0, 2.0];
    let (ls, us) = lu_nopivot(&a_s);
    let xs = back_sub(&us, &forward_sub(&ls, &b_s));
    println!("{:?}", xs.as_slice());
    println!("{:?}", (&b_s - &a_s * &xs).as_slice());      // residual
    let x_ok = a_s.clone().lu().solve(&b_s).unwrap();
    println!("{:?}", (&b_s - &a_s * &x_ok).as_slice());
}
