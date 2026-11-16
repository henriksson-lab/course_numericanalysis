// Q5: singular values of a matrix with two identical columns.
use nalgebra::DMatrix;

fn main() {
    let a = DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 2.0, 2.0, 3.0, 3.0]);
    println!("{:?}", a.singular_values().as_slice());
}
