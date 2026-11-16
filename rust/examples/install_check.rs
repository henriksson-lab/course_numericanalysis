use nalgebra::{matrix, vector};

fn main() {
    let a = matrix![2.0, 1.0;
                    1.0, 3.0];
    let b = vector![1.0, 2.0];
    let x = a.lu().solve(&b).unwrap();
    println!("{:?}", x.as_slice());   // should print [0.2, 0.6] (to rounding)
    println!("{}", 0.1 + 0.2 == 0.3); // false: see How computers store numbers
}
