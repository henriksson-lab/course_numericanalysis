// C2: explicit Euler vs backward Euler (with Newton, nalgebra LU) on the enzyme model.
use nalgebra::{DMatrix, DVector};

const KF: f64 = 100.0; const KR: f64 = 100.0; const KCAT: f64 = 1.0;

fn f(y: &DVector<f64>) -> DVector<f64> {                // y = (S, E, C, P)
    let b = KF * y[1] * y[0] - KR * y[2];
    DVector::from_vec(vec![-b, -b + KCAT * y[2], b - KCAT * y[2], KCAT * y[2]])
}

fn jac(y: &DVector<f64>) -> DMatrix<f64> {
    DMatrix::from_row_slice(4, 4, &[
        -KF * y[1], -KF * y[0],  KR,         0.0,
        -KF * y[1], -KF * y[0],  KR + KCAT,  0.0,
         KF * y[1],  KF * y[0], -KR - KCAT,  0.0,
         0.0,        0.0,        KCAT,       0.0])
}

fn euler(y0: &DVector<f64>, h: f64, n: usize) -> DVector<f64> {
    let mut y = y0.clone();
    for _ in 0..n { y = &y + h * f(&y); }
    y
}

// returns the last value and the total number of Newton iterations
fn beuler(y0: &DVector<f64>, h: f64, n: usize, tol: f64) -> (DVector<f64>, usize) {
    let mut y = y0.clone();
    let id = DMatrix::<f64>::identity(y.len(), y.len());
    let mut newton = 0;
    for _ in 0..n {
        let yold = y.clone();                            // Newton starts at y_n
        for _ in 0..20 {
            let g = &y - &yold - h * f(&y);
            let delta = (&id - h * jac(&y)).lu().solve(&(-g)).expect("singular I - hJ");
            y += &delta; newton += 1;
            if delta.norm() < tol { break; }
        }
    }
    (y, newton)
}

fn main() {
    let y0 = DVector::from_vec(vec![10.0, 1.0, 0.0, 0.0]);
    for h in [1e-3, 3e-3] {
        let y = euler(&y0, h, (10.0 / h) as usize);      // t = 10 s (3333 steps for 3 ms)
        println!("{}", y.iter().map(|v| format!("{:.6}", v)).collect::<Vec<_>>().join(" "));
    }
    for h in [1.0f64, 0.1, 0.01] {
        let n = (10.0 / h).round() as usize;
        let (y, newton) = beuler(&y0, h, n, 1e-10);
        println!("h = {:4.2}  P(10) = {:.6}  Newton its/step = {:.2}  E+C = {:.10}  S+C+P = {:.10}",
                 h, y[3], newton as f64 / n as f64, y[1] + y[2], y[0] + y[2] + y[3]);
    }
}
