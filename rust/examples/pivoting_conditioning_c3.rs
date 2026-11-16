// Vector and matrix norms with nalgebra, and the 2-norm as the largest stretch.
use nalgebra::{matrix, vector, Vector2};

fn main() {
    let v = vector![3.0, -4.0, 12.0];
    println!("{} {} {}", v.lp_norm(1), v.norm(), v.amax());
    let dx = vector![0.2, -0.2];
    let x = vector![1.0, 1.0];
    // relative error of x-hat = (1.2, 0.8) in the 1-, 2- and inf-norm
    println!("{} {} {}", dx.lp_norm(1) / x.lp_norm(1), dx.norm() / x.norm(), dx.amax() / x.amax());

    let a = matrix![1.0,  0.95;
                    0.95, 1.0];
    let abs_a = a.abs();
    // row_sum() = column sums -> 1-norm; column_sum() = row sums -> inf-norm
    println!("{} {} {}", abs_a.row_sum().max(), abs_a.column_sum().max(),
             a.singular_values().max());

    // stretch factors |A v| for 1001 unit vectors v on the circle
    let stretch: Vec<f64> = (0..=1000)
        .map(|i| {
            let th = 2.0 * std::f64::consts::PI * i as f64 / 1000.0;
            (a * Vector2::new(th.cos(), th.sin())).norm()
        })
        .collect();
    println!("{}", stretch.iter().cloned().fold(f64::MIN, f64::max)); // largest
    println!("{}", stretch.iter().cloned().fold(f64::MAX, f64::min)); // smallest
}

// Output:
// 19 13 12
// 0.2 0.2 0.2
// 1.95 1.95 1.95
// 1.95
// 0.05000000000000006
