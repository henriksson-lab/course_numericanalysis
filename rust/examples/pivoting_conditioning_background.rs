// Solve, norms and condition number with nalgebra
// (R: solve, norm, kappa; numpy: linalg.solve, linalg.norm, linalg.cond).
use nalgebra::{matrix, vector};

fn main() {
    let a = matrix![1.0,  0.95;
                    0.95, 1.0];
    let x = a.lu().solve(&vector![1.96, 1.94]).unwrap(); // LU with partial pivoting
    println!("{}  {}", x[0], x[1]);
    let v = vector![3.0, -4.0, 12.0];
    println!("{} {} {}", v.lp_norm(1), v.norm(), v.amax()); // 1-, 2- and inf-norm
    // nalgebra has no induced 1- and inf-norm. row_sum() adds up the rows, giving
    // the column sums (1-norm); column_sum() gives the row sums (inf-norm).
    let abs_a = a.abs();
    let s = a.singular_values();
    println!("{} {} {}", abs_a.row_sum().max(), abs_a.column_sum().max(), s.max());
    println!("{}", s.max() / s.min());                   // condition number (2-norm)
}

// Output:
// 1.1999999999999995  0.8000000000000005
// 19 13 12
// 1.95 1.95 1.95
// 38.999999999999964
