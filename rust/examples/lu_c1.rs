// C1: LU without pivoting, forward and back substitution, written out.
use nalgebra::{dmatrix, dvector, DMatrix, DVector};

fn lu_nopivot(a: &DMatrix<f64>) -> (DMatrix<f64>, DMatrix<f64>) {
    let n = a.nrows();
    let mut l = DMatrix::<f64>::identity(n, n);
    let mut u = a.clone();
    for k in 0..n - 1 {                           // pivot column k
        for i in k + 1..n {                       // rows below the pivot
            l[(i, k)] = u[(i, k)] / u[(k, k)];    // multiplier l_ik
            for j in k..n {                       // row i - l_ik * row k
                u[(i, j)] -= l[(i, k)] * u[(k, j)];
            }
        }
    }
    (l, u)
}

fn forward_sub(l: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let n = b.len();
    let mut y = DVector::<f64>::zeros(n);
    for i in 0..n {                               // top to bottom
        let s: f64 = (0..i).map(|j| l[(i, j)] * y[j]).sum();
        y[i] = (b[i] - s) / l[(i, i)];
    }
    y
}

fn back_sub(u: &DMatrix<f64>, y: &DVector<f64>) -> DVector<f64> {
    let n = y.len();
    let mut x = DVector::<f64>::zeros(n);
    for i in (0..n).rev() {                       // bottom to top
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
    let a = dmatrix![1.0, 1.0, 1.0; -13.0, -21.0, -27.0; 6.0, 18.0, 2.0];
    let (l, u) = lu_nopivot(&a);
    show(&l);
    show(&u);
    println!("{:e}", (&l * &u - &a).abs().max());
    for b in [dvector![1.0, -18.4, 11.6], dvector![1.0, -22.4, 7.6]] {   // two mussels
        let y = forward_sub(&l, &b);
        println!("{:?}", y.as_slice());
        println!("{:?}", back_sub(&u, &y).as_slice());
    }
    let a2 = dmatrix![2.0, 1.0, 1.0; 4.0, 3.0, 3.0; 8.0, 7.0, 9.0];
    let (l2, u2) = lu_nopivot(&a2);
    show(&l2);
    show(&u2);
    println!("{}", u2.diagonal().product());
    println!("{}", a2.clone().lu().determinant());
}
