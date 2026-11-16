// C1: qPCR standard curve by the normal equations; efficiency and unknowns.
use nalgebra::{DMatrix, DVector};

// scientific notation as R and Python print it: 2.2e+01 rather than Rust's 2.2e1
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let log_n: Vec<f64> = (1..=7).flat_map(|k| [k as f64; 3]).collect();
    let cq = DVector::from_vec(vec![
        34.52, 34.69, 34.48, 31.62, 31.37, 31.14, 28.05, 28.10, 28.07,
        24.54, 24.90, 24.68, 21.13, 20.81, 21.47, 17.89, 17.90, 18.09,
        14.71, 14.67, 14.73]);
    let a = DMatrix::from_fn(log_n.len(), 2, |i, j| if j == 0 { 1.0 } else { log_n[i] });
    let x = (a.transpose() * &a).lu().solve(&(a.transpose() * &cq)).unwrap();
    println!("{:.6} {:.6}", x[0], x[1]);
    let r = &cq - &a * &x;
    let atr = a.transpose() * &r;                      // A^T r, should be ~ 0
    println!("{} {}", sci(atr[0], 6), sci(atr[1], 6));
    println!("{:.7}", r.norm_squared());
    println!("{:.7}", 10f64.powf(-1.0 / x[1]) - 1.0);  // efficiency
    for cq_unk in [24.0, 30.5] {
        let l = (cq_unk - x[0]) / x[1];
        println!("{:.6} {:.4}", l, 10f64.powf(l));
    }
}
