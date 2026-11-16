// C2: error per iteration and the measured convergence rate.
use nalgebra::{DMatrix, DVector};

fn main() {
    let l = DMatrix::from_row_slice(4, 4, &[
        0.0, 1.5, 4.0, 5.0,
        0.3, 0.0, 0.0, 0.0,
        0.0, 0.6, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
    ]);
    // "exact" answer: many power steps (nalgebra has no nonsymmetric eigenvectors)
    let mut v = DVector::from_row_slice(&[0.25, 0.25, 0.25, 0.25]);
    for _ in 0..500 {
        v = &l * &v;
        v /= v.sum();
    }
    let mut ev: Vec<f64> = l.complex_eigenvalues().iter().map(|z| z.norm()).collect();
    ev.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let ratio = ev[1] / ev[0];

    let kmax = 60;
    let mut err_x = vec![0.0; kmax];
    let mut x = DVector::from_row_slice(&[0.0, 1.0, 0.0, 0.0]);
    for k in 0..kmax {
        let y = &l * &x;
        x = &y / y.sum();
        err_x[k] = (&x - &v).amax();             // max |x - v|
    }
    // least-squares slope of log10(error) against k for k = 10..50
    let ks: Vec<f64> = (10..=50).map(|k| k as f64).collect();
    let es: Vec<f64> = (10..=50).map(|k| err_x[k - 1].log10()).collect();
    let (mk, me) = (ks.iter().sum::<f64>() / 41.0, es.iter().sum::<f64>() / 41.0);
    let slope = ks.iter().zip(&es).map(|(k, e)| (k - mk) * (e - me)).sum::<f64>()
        / ks.iter().map(|k| (k - mk).powi(2)).sum::<f64>();
    println!("{}", slope);                          // measured slope
    println!("{}", ratio.log10());                  // predicted slope
    println!("{}", 1e-10f64.ln() / ratio.ln());
}
