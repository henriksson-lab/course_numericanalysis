// C2: log det and quadratic form of a 2000 x 2000 LD matrix by Cholesky (nalgebra).
use nalgebra::{DMatrix, DVector};

fn main() {
    let (n, rho) = (2000usize, 0.9_f64);
    let s = DMatrix::from_fn(n, n, |i, j| rho.powi((i as i32 - j as i32).abs()));
    println!("{}", s.clone().lu().determinant());       // underflows
    let l = s.clone().cholesky().expect("not SPD").l();
    let logdet: f64 = 2.0 * l.diagonal().iter().map(|v| v.ln()).sum::<f64>();
    println!("{}", logdet);
    println!("{}", (n - 1) as f64 * (1.0 - rho * rho).ln()); // exact
    let y = DVector::from_element(n, 1.0);
    let z = l.solve_lower_triangular(&y).unwrap();      // z = L^{-1} y
    println!("{}", z.norm_squared());                    // y^T S^{-1} y
    println!("{}", y.dot(&s.lu().solve(&y).unwrap()));   // check with a general solve
    let nf = n as f64;
    println!("{}", (2.0 + (nf - 2.0) * (1.0 + rho * rho) - 2.0 * (nf - 1.0) * rho)
        / (1.0 - rho * rho));                            // exact
}
