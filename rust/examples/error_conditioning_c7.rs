// Amplification of a perturbation in b, compared with the condition number.
// Uses the nalgebra crate (R: solve, kappa; numpy: linalg.solve, linalg.cond).
use nalgebra::{Matrix2, Vector2};

// Amplification (|dx|/|x|) / (|db|/|b|); .norm() is the Euclidean norm.
fn amplification(a: &Matrix2<f64>, b: &Vector2<f64>, db: &Vector2<f64>) -> f64 {
    let lu = a.lu(); // LU factorisation, reused for both solves
    let x = lu.solve(b).unwrap(); // unwrap: stop if A is singular
    let dx = lu.solve(&(b + db)).unwrap() - x;
    (dx.norm() / x.norm()) / (db.norm() / b.norm())
}

// 2-norm condition number = largest / smallest singular value
fn cond(a: &Matrix2<f64>) -> f64 {
    let s = a.singular_values();
    s.max() / s.min()
}

fn main() {
    let a = Matrix2::new(1.0, 0.95,
                         0.95, 1.0); // entries given row by row
    let b = Vector2::new(1.95, 1.95);
    let db = Vector2::new(0.01, -0.01);
    println!("{}", amplification(&a, &b, &db));
    println!("{}", cond(&a));

    let a2 = Matrix2::new(1.0, 0.5,
                          0.5, 1.0);
    println!("{}", amplification(&a2, &b, &db));
    println!("{}", cond(&a2));
}
