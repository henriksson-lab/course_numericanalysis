// C3: too-long elements give negative concentrations; a two-layer gel with and
// without a node at the layer interface.

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
    let uniform = |n: usize| (0..=n).map(|i| l * i as f64 / n as f64).collect::<Vec<f64>>();
    let minval = |u: &[f64]| u.iter().cloned().fold(f64::INFINITY, f64::min);
    for n in [2, 4, 6, 8, 10] {
        let x = uniform(n);
        let u = fem1d(&x, &vec![dd; n], &vec![k; n], &vec![0.0; n], 1.0, false);
        let ul = fem1d(&x, &vec![dd; n], &vec![k; n], &vec![0.0; n], 1.0, true);
        println!("N = {:2}  h/lambda = {:4.1}  min u = {:9.2e}  lumped min = {:.2e}",
                 n, l / n as f64 / lam, minval(&u), minval(&ul));
    }
    // cell-free top layer (k = 0) of thickness c, cells below
    let c = 100.0;
    let uc = 1.0 / (1.0 + c / lam * ((l - c) / lam).tanh());
    println!("{}", uc);
    for n in [50, 49, 100, 99, 200, 199] {
        let x = uniform(n);
        let kk: Vec<f64> = (0..n).map(|e| if (x[e] + x[e + 1]) / 2.0 < c { 0.0 } else { k }).collect();
        let u = fem1d(&x, &vec![dd; n], &kk, &vec![0.0; n], 1.0, false);
        let e = x.iter().position(|&xi| xi >= c).unwrap();   // linear interpolation at x = c
        let u100 = if (x[e] - c).abs() < 1e-9 { u[e] }
                   else { u[e - 1] + (u[e] - u[e - 1]) * (c - x[e - 1]) / (x[e] - x[e - 1]) };
        let node = x.iter().any(|&xi| (xi - c).abs() < 1e-9);
        println!("N = {:3}  node at 100: {:5}  u(100) = {:.4}  error = {:9.2e}", n, node, u100, u100 - uc);
    }
}
