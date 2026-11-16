// Isotope mixing model: solve A f = b with nalgebra's LU (partial pivoting).
use nalgebra::{dmatrix, dvector};

fn main() {
    let a = dmatrix![  1.0,   1.0,   1.0;      // fractions sum to 1
                     -13.0, -21.0, -27.0;      // d13C: Spartina, phytoplankton, upland C3
                       6.0,  18.0,   2.0];     // d34S
    let b = dvector![1.0, -18.4, 11.6];        // one ribbed mussel
    let f = a.clone().lu();                    // factorise once ...
    let x = f.solve(&b).expect("singular");    // ... substitute as often as needed
    println!("{:?}", x.as_slice());
}
