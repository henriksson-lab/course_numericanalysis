// C2: kappa(A^T A) = kappa(A)^2, and the Vandermonde matrix of the BCA standards.
use nalgebra::DMatrix;

// scientific notation as R and Python print it: 2.2e+01 rather than Rust's 2.2e1
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn kappa(a: &DMatrix<f64>) -> f64 {
    let sv = a.clone().singular_values();
    sv.max() / sv.min()
}

// columns 1, x, x^2, ..., x^d
fn vander(x: &[f64], d: usize) -> DMatrix<f64> {
    DMatrix::from_fn(x.len(), d + 1, |i, j| x[i].powi(j as i32))
}

fn main() {
    let log_n: Vec<f64> = (1..=7).flat_map(|k| [k as f64; 3]).collect();
    let a = vander(&log_n, 1);
    println!("{} {} {}", kappa(&a), kappa(&(a.transpose() * &a)), kappa(&a).powi(2));
    let conc: Vec<f64> = [0., 25., 125., 250., 500., 750., 1000., 1500., 2000.]
        .iter().flat_map(|&c| [c; 2]).collect();
    let conc_mg: Vec<f64> = conc.iter().map(|c| c / 1000.0).collect();
    for d in 1..=4 {
        let k = kappa(&vander(&conc, d));
        println!("{} {} {} {:.3}", d, sci(k, 2), sci(k * k, 2), kappa(&vander(&conc_mg, d)));
    }
}
