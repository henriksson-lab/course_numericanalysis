// Q7: solve A x = b without forming the inverse.
use nalgebra::{dmatrix, dvector};

fn main() {
    let a = dmatrix![1.0, 1.0, 1.0; -26.0, -10.0, -18.0; 2.0, 6.0, 14.0];
    let b = dvector![1.0, -19.6, 5.6];
    let x = a.lu().solve(&b).expect("singular");   // not a.try_inverse().unwrap() * b
    println!("{:?}", x.as_slice());
}
