// Digits lost when unmixing 6 dyes with Gaussian spectra (nalgebra LU and SVD).
use nalgebra::{DMatrix, DVector};

// detector i sits at the peak of dye i; entry (i, j) = emission of dye j at detector i
fn panel(n: usize, spacing: f64, width: f64) -> DMatrix<f64> {
    let peak = |i: usize| 500.0 + spacing * i as f64;
    DMatrix::from_fn(n, n, |i, j| (-(peak(i) - peak(j)).powi(2) / (2.0 * width * width)).exp())
}

fn main() {
    let xt = DVector::from_vec(vec![5.0, 1.0, 3.0, 2.0, 4.0, 1.0]);
    for s in [60.0, 30.0, 20.0, 10.0, 5.0] {
        let a = panel(6, s, 30.0);
        let sv = a.singular_values();
        let k = sv.max() / sv.min();
        let xh = a.clone().lu().solve(&(&a * &xt)).unwrap();
        let err = (&xh - &xt).norm() / xt.norm();
        println!("spacing {:2} nm  kappa = {:.1e}  error = {:.1e}  digits = {:4.1}  16 - log10(kappa) = {:4.1}",
                 s, k, err, -err.log10(), 16.0 - k.log10());
    }
}

// Output:
// spacing 60 nm  kappa = 1.6e0  error = 6.8e-17  digits = 16.2  16 - log10(kappa) = 15.8
// spacing 30 nm  kappa = 2.9e1  error = 9.3e-16  digits = 15.0  16 - log10(kappa) = 14.5
// spacing 20 nm  kappa = 7.1e2  error = 3.8e-14  digits = 13.4  16 - log10(kappa) = 13.2
// spacing 10 nm  kappa = 6.4e5  error = 3.8e-11  digits = 10.4  16 - log10(kappa) = 10.2
// spacing  5 nm  kappa = 7.2e8  error = 7.9e-10  digits =  9.1  16 - log10(kappa) =  7.1
