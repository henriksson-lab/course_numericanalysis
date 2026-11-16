// Quiz Q8: the QR algorithm on the period-2 transition matrix, unshifted and shifted.
use nalgebra::Matrix2;

fn show(m: &Matrix2<f64>) {
    for i in 0..2 { println!("{:10.6} {:10.6}", m[(i, 0)], m[(i, 1)]); }
}

fn main() {
    let p: Matrix2<f64> = Matrix2::new(0.0, 1.0,
                                       1.0, 0.0);
    let qr = p.qr();
    show(&(qr.r() * qr.q()));                         // one unshifted step
    let mu = 0.5;
    let eye = Matrix2::<f64>::identity();
    let mut p1 = p;
    for _ in 0..10 {                                  // shifted steps, shift 0.5
        let qr = (p1 - eye * mu).qr();
        p1 = qr.r() * qr.q() + eye * mu;
    }
    show(&p1);
}
