// C3: 286 consumers against the same 4-source mixing matrix.
use nalgebra::{dmatrix, DMatrix};

fn main() {
    let s = dmatrix![  1.0,   1.0,   1.0,   1.0;
                     -13.0, -21.0, -27.0, -17.0;   // d13C: Spartina, phytoplankton, upland C3, microalgae
                       6.0,   8.0,   4.0,   5.0;   // d15N
                       6.0,  18.0,   2.0,  10.0];  // d34S
    // all diets with fractions in steps of 0.1 (286 consumers), one per column
    let mut cols: Vec<f64> = Vec::new();
    for i in 0..=10 {
        for j in 0..=(10 - i) {
            for k in 0..=(10 - i - j) {
                cols.extend([i, j, k, 10 - i - j - k].map(|v| v as f64 / 10.0));
            }
        }
    }
    let m = cols.len() / 4;
    let ftrue = DMatrix::from_column_slice(4, m, &cols);
    let mut b = &s * &ftrue;
    for c in 0..m {                                // deterministic "measurement noise"
        let t = (c + 1) as f64;
        b[(1, c)] += 0.2 * t.sin();
        b[(2, c)] += 0.2 * t.cos();
        b[(3, c)] += 0.2 * (2.0 * t).sin();
    }
    let x = s.lu().solve(&b).expect("singular");  // one factorisation, m solves
    println!("{m}");
    let colsum_err = x.column_iter().map(|c| (c.sum() - 1.0).abs()).fold(0.0, f64::max);
    println!("{:e}", colsum_err);
    let mins: Vec<f64> = x.column_iter().map(|c| c.min()).collect();
    println!("{}", mins.iter().filter(|&&v| v < 0.0).count());  // a negative fraction
    println!("{}", mins.iter().filter(|&&v| v < -0.1).count());
    println!("[{:.3}, {:.3}]", x.min(), x.max());
    println!("{}", (&x - &ftrue).abs().max());
}
