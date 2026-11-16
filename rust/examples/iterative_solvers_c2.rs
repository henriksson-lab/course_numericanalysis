// rust/examples/iterative_solvers_c2.rs  —  cargo run --example iterative_solvers_c2
// C2: spectral radius of the Jacobi and Gauss-Seidel iteration matrices (dense, nalgebra).
use nalgebra::DMatrix;
use std::f64::consts::PI;

fn lap2_dense(m: usize) -> DMatrix<f64> {
    let n = m * m;
    DMatrix::from_fn(n, n, |k, l| {
        let (i, j, p, q) = (k % m, k / m, l % m, l / m);
        if k == l { 4.0 }
        else if (i == p && j.abs_diff(q) == 1) || (j == q && i.abs_diff(p) == 1) { -1.0 }
        else { 0.0 }
    })
}
fn rho(t: &DMatrix<f64>) -> f64 {                   // largest |eigenvalue|
    t.complex_eigenvalues().iter().map(|z| z.norm()).fold(0.0, f64::max)
}

fn main() {
    let m = 7;
    let a = lap2_dense(m);
    let id = DMatrix::<f64>::identity(m * m, m * m);
    let tj = &id - &a / 4.0;                        // I - D^{-1} A, D = 4 I
    let low = a.lower_triangle();
    let tgs = &id - low.solve_lower_triangular(&a).unwrap();
    let (rj, rg) = (rho(&tj), rho(&tgs));
    println!("{} {} {} {}", rj, (PI / 8.0).cos(), rg, (PI / 8.0).cos().powi(2));
    println!("{} {}", 10f64.ln() / -rj.ln(), 10f64.ln() / -rg.ln());   // iterations per digit
    let b2 = DMatrix::from_row_slice(2, 2, &[1.0, 2.0, 2.0, 1.0]);
    println!("{}", rho(&(DMatrix::<f64>::identity(2, 2) - b2)));     // the 2 x 2 failure
}
