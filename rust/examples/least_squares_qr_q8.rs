// Q8: an intercept, a dose column and a column made from both.
use nalgebra::DMatrix;

fn main() {
    let dose = [0.0, 1.0, 2.0, 4.0, 8.0];
    let a = DMatrix::from_fn(5, 3, |i, j| match j {
        0 => 1.0,
        1 => dose[i],
        _ => dose[i] / 2.0 + 1.0,
    });
    let sv = a.singular_values();
    println!("{}", a.rank(1e-7 * sv.max()));
}
