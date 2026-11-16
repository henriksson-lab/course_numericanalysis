// Background: the toggle switch F, its Jacobian J, and one Newton step.
use nalgebra::{Matrix2, Vector2};

const ALPHA: f64 = 3.0;
const N: f64 = 2.0;

// x = (u, v); returns (F1, F2)
fn f(x: &Vector2<f64>) -> Vector2<f64> {
    let (u, v) = (x[0], x[1]);
    Vector2::new(ALPHA / (1.0 + v.powf(N)) - u,
                 ALPHA / (1.0 + u.powf(N)) - v)
}

// row i = equation i, column j = variable j
fn jac(x: &Vector2<f64>) -> Matrix2<f64> {
    let (u, v) = (x[0], x[1]);
    Matrix2::new(-1.0, -ALPHA * N * v.powf(N - 1.0) / (1.0 + v.powf(N)).powi(2),
                 -ALPHA * N * u.powf(N - 1.0) / (1.0 + u.powf(N)).powi(2), -1.0)
}

fn main() {
    let x0 = Vector2::new(3.0, 0.5);
    println!("{:?}", f(&x0).as_slice());
    println!("{}", jac(&x0));            // nalgebra prints matrices as a box
    // one Newton step = one linear solve (LU); solve() returns None if J is singular
    let delta = jac(&x0).lu().solve(&(-f(&x0))).expect("singular Jacobian");
    println!("{:?}", (x0 + delta).as_slice());
}
