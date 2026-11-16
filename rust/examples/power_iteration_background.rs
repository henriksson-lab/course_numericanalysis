// Background: eigenvalues of the fish Leslie matrix with nalgebra.
// nalgebra has no general (nonsymmetric) eigenvector routine, so we only get
// the eigenvalues here; the eigenvector comes from power iteration (C1).
use nalgebra::DMatrix;

fn main() {
    let l = DMatrix::from_row_slice(4, 4, &[
        0.0, 1.5, 4.0, 5.0,      // offspring per fish of each age
        0.3, 0.0, 0.0, 0.0,      // survival to the next age
        0.0, 0.6, 0.0, 0.0,
        0.0, 0.0, 0.5, 0.0,
    ]);
    let mut ev: Vec<_> = l.complex_eigenvalues().iter().cloned().collect();
    ev.sort_by(|a, b| b.norm().partial_cmp(&a.norm()).unwrap()); // largest |lambda| first
    for z in &ev {
        println!("{:.7} {:+.7}i   |lambda| = {:.7}", z.re, z.im, z.norm());
    }
}
