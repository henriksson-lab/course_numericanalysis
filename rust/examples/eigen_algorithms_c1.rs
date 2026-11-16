// C1: the unshifted QR algorithm on a symmetric 3 x 3 matrix.
use nalgebra::Matrix3;

// scientific notation as R and Python print it: 1.33e-02 rather than Rust's 1.33e-2
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let a = Matrix3::new(4.0, 1.0, 0.0,
                         1.0, 3.0, 1.0,
                         0.0, 1.0, 2.0);
    let mut lam: Vec<f64> = a.symmetric_eigen().eigenvalues.iter().cloned().collect();
    lam.sort_by(|x, y| y.partial_cmp(x).unwrap());
    let mut ak = a;
    let mut qk = Matrix3::<f64>::identity();
    for k in 1..=40 {
        let qr = ak.qr();
        qk *= qr.q();                       // accumulate Q1 Q2 ... Qk
        ak = qr.r() * qr.q();
        if k % 10 == 0 {
            println!("k = {:2}  |a21| = {}  ratio^k = {}  |a32| = {}  ratio^k = {}", k,
                     sci(ak[(1, 0)].abs(), 2), sci((lam[1] / lam[0]).powi(k), 2),
                     sci(ak[(2, 1)].abs(), 2), sci((lam[2] / lam[1]).powi(k), 2));
        }
    }
    let d: Vec<String> = (0..3).map(|i| sci(ak[(i, i)] - lam[i], 6)).collect();
    println!("{}", d.join(" "));            // eigenvalues on the diagonal
    println!("{}", sci((qk.transpose() * a * qk - ak).amax(), 6)); // A_k = Q^T A Q
    // simultaneous iteration: Q1...Q5 against the Q of A^5 (up to signs)
    let (mut ak, mut qk) = (a, Matrix3::<f64>::identity());
    for _ in 0..5 {
        let qr = ak.qr();
        qk *= qr.q();
        ak = qr.r() * qr.q();
    }
    let q5 = (a * a * a * a * a).qr().q();
    println!("{}", sci((qk.abs() - q5.abs()).amax(), 6));
}
