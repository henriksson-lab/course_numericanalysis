// C2: shifted QR algorithm with deflation; no shift, Rayleigh shift, Wilkinson shift.
use nalgebra::DMatrix;

fn shift_of(a: &DMatrix<f64>, kind: &str) -> f64 {
    let n = a.nrows();
    match kind {
        "none" => 0.0,
        "rayleigh" => a[(n - 1, n - 1)],
        _ => {                                            // Wilkinson
            let (p, b, c) = (a[(n - 2, n - 2)], a[(n - 1, n - 2)], a[(n - 1, n - 1)]);
            let d = (p - c) / 2.0;
            let s = if d >= 0.0 { 1.0 } else { -1.0 };
            c - s * b * b / (d.abs() + (d * d + b * b).sqrt())
        }
    }
}

fn qr_eig(mut a: DMatrix<f64>, kind: &str, tol: f64, maxit: usize) -> (Vec<f64>, usize) {
    let mut lam: Vec<f64> = Vec::new();
    let mut steps = 0;
    while a.nrows() > 1 && steps < maxit {
        let n = a.nrows();
        let mu = shift_of(&a, kind);
        let eye = DMatrix::<f64>::identity(n, n);
        let qr = (&a - &eye * mu).qr();
        a = qr.r() * qr.q() + eye * mu;
        steps += 1;
        if a[(n - 1, n - 2)].abs() < tol {               // deflate
            lam.insert(0, a[(n - 1, n - 1)]);
            a = a.view((0, 0), (n - 1, n - 1)).into_owned();
        }
    }
    let mut vals: Vec<f64> = a.diagonal().iter().cloned().collect();
    vals.extend(lam);
    (vals, steps)
}

fn main() {
    let a = DMatrix::from_row_slice(3, 3, &[4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0]);
    let b = DMatrix::from_row_slice(2, 2, &[2.0, 1.0, 1.0, 2.0]);
    for m in [&a, &b] {
        for kind in ["none", "rayleigh", "wilkinson"] {
            let (vals, steps) = qr_eig(m.clone(), kind, 1e-12, 500);
            let v: Vec<String> = vals.iter().map(|x| format!("{:.10}", x)).collect();
            println!("{:<9} steps = {:3}  values = {}", kind, steps, v.join(" "));
        }
    }
}
