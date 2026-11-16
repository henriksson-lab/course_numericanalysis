// Q11: no stiff library solver in Rust, so check the hand-written backward Euler
// (log-spaced steps from 1e-6 to 40) by refining the grid and by the conservation law.
use nalgebra::{DMatrix, DVector};

fn f(y: &DVector<f64>) -> DVector<f64> {
    let (r1, r2, r3) = (0.04 * y[0], 3e7 * y[1] * y[1], 1e4 * y[1] * y[2]);
    DVector::from_vec(vec![-r1 + r3, r1 - r2 - r3, r2])
}
fn jac(y: &DVector<f64>) -> DMatrix<f64> {
    DMatrix::from_row_slice(3, 3, &[
        -0.04, 1e4 * y[2],                 1e4 * y[1],
         0.04, -6e7 * y[1] - 1e4 * y[2], -1e4 * y[1],
         0.0,   6e7 * y[1],                0.0])
}

fn main() {
    let id = DMatrix::<f64>::identity(3, 3);
    for m in [100, 1000, 10000] {                       // steps per decade
        let n = (m as f64 * (40f64.log10() + 6.0)).round() as usize;
        let times: Vec<f64> = std::iter::once(0.0)
            .chain((0..=n).map(|k| 1e-6 * (40e6f64).powf(k as f64 / n as f64))).collect();
        let mut y = DVector::from_vec(vec![1.0, 0.0, 0.0]);
        for k in 1..times.len() {
            let (h, yold) = (times[k] - times[k - 1], y.clone());
            for _ in 0..20 {
                let delta = (&id - h * jac(&y)).lu().solve(&(-(&y - &yold - h * f(&y)))).unwrap();
                y += &delta;
                if delta.norm() < 1e-14 { break; }
            }
        }
        println!("{:5} steps per decade: y2(40) = {:.7e}  sum - 1 = {:.1e}", m, y[1], y.sum() - 1.0);
    }
}
