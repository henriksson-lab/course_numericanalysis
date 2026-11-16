// Background: eigenvalues of a symmetric and a nonsymmetric matrix, and one QR step.
use nalgebra::{Matrix2, Matrix3};

fn main() {
    let a = Matrix3::new(4.0, 1.0, 0.0,
                         1.0, 3.0, 1.0,
                         0.0, 1.0, 2.0);                 // symmetric 3 x 3
    let mut ev: Vec<f64> = a.symmetric_eigen().eigenvalues.iter().cloned().collect();
    ev.sort_by(|x, y| y.partial_cmp(x).unwrap());         // unsorted: sort decreasing
    println!("{:.6?}", ev);
    let j = Matrix2::new(-0.1, -1.0,
                          0.5,  0.0);                     // predator-prey Jacobian
    for z in j.complex_eigenvalues().iter() {
        println!("{:.7} {:+.7}i", z.re, z.im);
    }
    let qr = a.qr();                                      // one QR step: A1 = R Q
    let a1 = qr.r() * qr.q();
    for i in 0..3 {
        println!("{:8.4} {:8.4} {:8.4}", a1[(i, 0)], a1[(i, 1)], a1[(i, 2)]);
    }
}
