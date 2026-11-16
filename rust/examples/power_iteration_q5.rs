// Q5: power iteration on the biennial plant after 100 and 101 steps.
use nalgebra::{DMatrix, DVector};

fn main() {
    let l = DMatrix::from_row_slice(2, 2, &[0.0, 4.0, 0.25, 0.0]);
    for steps in [100, 101] {
        let mut x = DVector::from_row_slice(&[1.0, 0.0]);
        for _ in 0..steps {
            let y = &l * &x;
            x = &y / y.sum();
        }
        println!("{} {:?}", steps, x.as_slice());
    }
}
