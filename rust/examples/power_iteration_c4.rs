// C4: stationary distribution of the cell-state chain three ways; the slow chain.
use nalgebra::{DMatrix, RowDVector};

fn main() {
    let p = DMatrix::from_row_slice(3, 3, &[
        0.80, 0.15, 0.05,
        0.05, 0.90, 0.05,
        0.05, 0.10, 0.85,
    ]);
    let mut d = RowDVector::from_row_slice(&[1.0, 0.0, 0.0]);   // row vector
    for _ in 0..100 {
        d = &d * &p;
    }
    println!("{:?}", d.as_slice());
    // linear system: (P^T - I) pi = 0 with the last equation replaced by sum(pi) = 1
    let mut m = p.transpose() - DMatrix::<f64>::identity(3, 3);
    m.row_mut(2).fill(1.0);
    let pi = m.lu().solve(&nalgebra::DVector::from_row_slice(&[0.0, 0.0, 1.0])).unwrap();
    println!("{:?}", pi.as_slice());

    let pslow = DMatrix::<f64>::identity(3, 3) + 0.01 * (&p - DMatrix::<f64>::identity(3, 3));
    let mut d = RowDVector::from_row_slice(&[1.0, 0.0, 0.0]);
    let mut change: Vec<f64> = Vec::new();
    let mut k = 0;
    while k < 100_000 {
        k += 1;
        let d_new = &d * &pslow;
        change.push((&d_new - &d).amax());
        d = d_new;
        if *change.last().unwrap() < 1e-8 {
            break;
        }
    }
    println!("{}", k);
    let n = change.len();
    let c_est = change[n - 1] / change[n - 2];
    println!("{}", c_est);
    println!("{:e}", change[n - 1] * c_est / (1.0 - c_est));        // estimated error
    let truth = RowDVector::from_row_slice(&[0.2, 0.55, 0.25]);
    println!("{:e}", (&d - truth).amax());                           // true error
}
