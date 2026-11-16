// rust/examples/phase_plane_c2.rs  —  cargo run --example phase_plane_c2
// C2: classify fixed points from the eigenvalues of a 2x2 Jacobian (nalgebra).
use nalgebra::{Complex, Matrix2, Vector2};

const ALPHA: f64 = 3.0;
const N: f64 = 2.0;

fn f(x: &Vector2<f64>) -> Vector2<f64> {
    let (u, v) = (x[0], x[1]);
    Vector2::new(ALPHA / (1.0 + v.powf(N)) - u, ALPHA / (1.0 + u.powf(N)) - v)
}

fn jac(x: &Vector2<f64>) -> Matrix2<f64> {
    let (u, v) = (x[0], x[1]);
    Matrix2::new(-1.0, -ALPHA * N * v.powf(N - 1.0) / (1.0 + v.powf(N)).powi(2),
                 -ALPHA * N * u.powf(N - 1.0) / (1.0 + u.powf(N)).powi(2), -1.0)
}

// eigenvalues sorted by real part, then imaginary part (largest first)
fn eig(a: &Matrix2<f64>) -> Vec<Complex<f64>> {
    let mut lam: Vec<Complex<f64>> = a.complex_eigenvalues().iter().cloned().collect();
    lam.sort_by(|p, q| (-p.re, -p.im).partial_cmp(&(-q.re, -q.im)).unwrap());
    lam
}

fn classify(a: &Matrix2<f64>, tol: f64) -> String {
    let lam = eig(a);
    if lam.iter().any(|z| z.re.abs() < tol) { return "undecided (Re = 0)".into(); }
    if lam.iter().any(|z| z.re > 0.0) && lam.iter().any(|z| z.re < 0.0) { return "saddle".into(); }
    let st = if lam.iter().all(|z| z.re < 0.0) { "stable" } else { "unstable" };
    let ty = if lam.iter().any(|z| z.im.abs() > tol) { "spiral" } else { "node" };
    format!("{} {}", st, ty)
}

fn fmt(lam: &[Complex<f64>]) -> String {
    lam.iter().map(|z| format!("{:8.4}{:+.4}i", z.re, z.im + 0.0)).collect::<Vec<_>>().join(" ")
}

fn fd_jacobian(f: impl Fn(&Vector2<f64>) -> Vector2<f64>, x: &Vector2<f64>) -> Matrix2<f64> {
    let h = f64::EPSILON.sqrt();
    let fx = f(x);
    let mut jm = Matrix2::zeros();
    for j in 0..2 {
        let hj = h * x[j].abs().max(1.0);
        let mut xh = *x;
        xh[j] += hj;
        jm.set_column(j, &((f(&xh) - fx) / hj));
    }
    jm
}

fn main() {
    let states = [("U on", Vector2::new(2.618034, 0.381966)),
                  ("middle", Vector2::new(1.213412, 1.213412)),
                  ("V on", Vector2::new(0.381966, 2.618034))];
    let mut mats = vec![("A1", Matrix2::new(-3.0, 1.0, 1.0, -3.0)),
                        ("A2", Matrix2::new(1.0, 2.0, 2.0, 1.0)),
                        ("A3", Matrix2::new(-1.0, -4.0, 1.0, -1.0)),
                        ("A4", Matrix2::new(0.0, 1.0, -4.0, 0.0))];
    for (name, x) in &states { mats.push((name, jac(x))); }
    for (name, a) in &mats {
        println!("{:<7} {}   {}", name, fmt(&eig(a)), classify(a, 1e-8));
    }
    for (name, x) in &states {
        let (e1, e2) = (eig(&fd_jacobian(f, x)), eig(&jac(x)));
        let d = e1.iter().zip(&e2).map(|(p, q)| (p - q).norm()).fold(0.0, f64::max);
        println!("{:<7} eigenvalues from FD Jacobian differ by {:.1e}", name, d);
    }
}
