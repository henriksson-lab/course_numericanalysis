// Background: steady advection-diffusion u c' = D c'' on [0, L], c(0) = 0, c(L) = 1,
// with central and upwind differences for c'. The tridiagonal system is solved with the
// Thomas algorithm (R and Python use a dense solve; same numbers).
const L: f64 = 500.0; // um
const U: f64 = 100.0; // um/s, Pe = u L / D = 50

fn exact(x: f64, d: f64) -> f64 {
    (U * x / d).exp_m1() / (U * L / d).exp_m1()
}

// a = sub-, b = main, c = super-diagonal, d = right-hand side
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec());
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    x
}

fn steady(n: usize, central: bool, d: f64) -> Vec<f64> {
    let dx = L / n as f64;
    let p = U * dx / d; // cell Peclet number
    let m = n - 1;
    let (lo, di, up) = if central { (1.0 + p / 2.0, -2.0, 1.0 - p / 2.0) } else { (1.0 + p, -2.0 - p, 1.0) };
    let mut rhs = vec![0.0; m];
    rhs[m - 1] = -up; // c(L) = 1 moves to the right side
    let inner = thomas(&vec![lo; m], &vec![di; m], &vec![up; m], &rhs);
    let mut c = vec![0.0];
    c.extend(inner);
    c.push(1.0);
    c
}

fn show(v: &[f64]) -> String {
    v.iter().map(|x| format!("{:6.3}", x)).collect::<Vec<_>>().join(" ")
}

fn main() {
    let x: Vec<f64> = (0..=10).map(|i| i as f64 * L / 10.0).collect();
    println!("{}", show(&steady(10, true, 1000.0)));
    println!("{}", show(&steady(10, false, 1000.0)));
    println!("{}", show(&x.iter().map(|&xi| exact(xi, 1000.0)).collect::<Vec<_>>()));
}
