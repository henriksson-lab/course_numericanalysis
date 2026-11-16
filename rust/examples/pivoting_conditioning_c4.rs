// Condition number of the two-dye matrix [1 c; c 1] from nalgebra's singular values.
use nalgebra::matrix;

fn main() {
    for c in [0.5, 0.9, 0.95, 0.99, 0.999] {
        let a = matrix![1.0, c;
                        c,   1.0];
        let s = a.singular_values();
        println!("c = {:5.3}  formula = {:7.1}  cond = {:7.1}",
                 c, (1.0 + c) / (1.0 - c), s.max() / s.min());
    }
}

// Output:
// c = 0.500  formula =     3.0  cond =     3.0
// c = 0.900  formula =    19.0  cond =    19.0
// c = 0.950  formula =    39.0  cond =    39.0
// c = 0.990  formula =   199.0  cond =   199.0
// c = 0.999  formula =  1999.0  cond =  1999.0
