// Background: a 1D finite element solver for -(D u')' + k u = f with linear
// elements on any mesh, hand-written (assembly loop + Thomas algorithm).
// Model: drug diffusing into a 1 mm hydrogel slab with uptake by cells.

// Thomas algorithm (Structure: symmetric, banded, sparse); a[0] and c[n-1] unused
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

/// Assemble the tridiagonal matrix (diagonal, off-diagonal) and load vector.
/// d, k, f hold one value per element.
fn assemble(x: &[f64], d: &[f64], k: &[f64], f: &[f64], lumped: bool) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = x.len();
    let (wd, wo) = if lumped { (0.5, 0.0) } else { (1.0 / 3.0, 1.0 / 6.0) };
    let (mut di, mut off, mut b) = (vec![0.0; n], vec![0.0; n - 1], vec![0.0; n]);
    for e in 0..n - 1 {                         // element e: nodes e and e + 1
        let h = x[e + 1] - x[e];
        let kd = d[e] / h + k[e] * h * wd;
        let ko = -d[e] / h + k[e] * h * wo;
        di[e] += kd;
        di[e + 1] += kd;
        off[e] += ko;
        b[e] += f[e] * h / 2.0;
        b[e + 1] += f[e] * h / 2.0;
    }
    (di, off, b)
}

/// Solve with u(x[0]) = u0 and no flux at the last node.
fn fem1d(x: &[f64], d: &[f64], k: &[f64], f: &[f64], u0: f64, lumped: bool) -> Vec<f64> {
    let (di, off, b) = assemble(x, d, k, f, lumped);
    let mut rhs = b[1..].to_vec();
    rhs[0] -= off[0] * u0;                     // move the known u0 to the right
    let mut c = off[1..].to_vec();
    c.push(0.0);
    let mut u = vec![u0];
    u.extend(thomas(&off, &di[1..], &c, &rhs));
    u
}

fn main() {
    let (l, dd, k) = (1000.0, 100.0, 0.04);
    let lam = f64::sqrt(dd / k);
    let exact = |x: f64| ((l - x) / lam).cosh() / (l / lam).cosh();
    let n = 10;                                // elements
    let x: Vec<f64> = (0..=n).map(|i| l * i as f64 / n as f64).collect();
    let u = fem1d(&x, &vec![dd; n], &vec![k; n], &vec![0.0; n], 1.0, false);
    println!("    x    fem  exact");
    for i in 0..5 {
        println!("{:5.0} {:6.4} {:6.4}", x[i], u[i], exact(x[i]));
    }
    let err = x.iter().zip(&u).map(|(xi, ui)| (ui - exact(*xi)).abs()).fold(0.0, f64::max);
    println!("{}", err);
}
