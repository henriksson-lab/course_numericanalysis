// C3: complex eigenvalues of a real matrix; the real QR algorithm cannot triangularise it.
use nalgebra::Matrix2;

fn main() {
    let j: Matrix2<f64> = Matrix2::new(-0.1, -1.0,
                                    0.5,  0.0);
    let ev = j.complex_eigenvalues();
    for z in ev.iter() {
        println!("{:.7} {:+.7}i", z.re, z.im);
    }
    println!("{:.6}", 2.0 * std::f64::consts::PI / ev[0].im.abs()); // period
    let mut jk = j;
    for k in 1..=200 {
        let qr = jk.qr();
        jk = qr.r() * qr.q();
        if [50, 51, 200].contains(&k) {
            println!("k = {:3}  a21 = {:8.4}  diag = {:8.4} {:8.4}",
                     k, jk[(1, 0)], jk[(0, 0)], jk[(1, 1)]);
        }
    }
    let (tr, dt) = (jk.trace(), jk.determinant());  // eigenvalues of the 2 x 2 block
    let disc = tr * tr - 4.0 * dt;                   // negative: complex pair
    println!("{:.7} {:+.7}i", tr / 2.0, (-disc).sqrt() / 2.0);
}
