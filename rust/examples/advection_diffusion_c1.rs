// C1: wiggles and smearing in the steady problem: max error vs grid, a wiggle check,
// observed orders, and the identity upwind(D) = central(D + u dx / 2).
const L: f64 = 500.0;
const U: f64 = 100.0;
const D: f64 = 1000.0;

fn exact(x: f64) -> f64 {
    (U * x / D).exp_m1() / (U * L / D).exp_m1()
}

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
    let p = U * (L / n as f64) / d;
    let m = n - 1;
    let (lo, di, up) = if central { (1.0 + p / 2.0, -2.0, 1.0 - p / 2.0) } else { (1.0 + p, -2.0 - p, 1.0) };
    let mut rhs = vec![0.0; m];
    rhs[m - 1] = -up;
    let mut c = vec![0.0];
    c.extend(thomas(&vec![lo; m], &vec![di; m], &vec![up; m], &rhs));
    c.push(1.0);
    c
}

fn wiggly(c: &[f64]) -> bool {
    c.windows(2).any(|p| p[1] - p[0] < 0.0) // exact profile only increases
}

fn sci(x: f64, digits: usize) -> String { // 1.60e-01 as R/Python print it
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let ns = [10, 20, 25, 40, 50, 100, 200, 400, 800];
    let (mut ec, mut eu) = (vec![], vec![]);
    for &n in &ns {
        let x: Vec<f64> = (0..=n).map(|i| i as f64 * L / n as f64).collect();
        let (cc, cu) = (steady(n, true, D), steady(n, false, D));
        let err = |c: &[f64]| c.iter().zip(&x).map(|(ci, xi)| (ci - exact(*xi)).abs()).fold(0.0, f64::max);
        ec.push(err(&cc));
        eu.push(err(&cu));
        println!("{:4} {:6.3} {:5.2} {:>9} {:<5} {:>9} {:<5}", n, L / n as f64, U * L / n as f64 / D,
                 sci(*ec.last().unwrap(), 2), wiggly(&cc).to_string().to_uppercase(),
                 sci(*eu.last().unwrap(), 2), wiggly(&cu).to_string().to_uppercase());
    }
    let k = ns.len() - 1;
    println!("order: central {:.2}, upwind {:.2}", (ec[k - 1] / ec[k]).log2(), (eu[k - 1] / eu[k]).log2());
    for n in [10, 50, 200] {
        let dx = L / n as f64;
        let a = steady(n, false, D);
        let b = steady(n, true, D + U * dx / 2.0);
        let diff = a.iter().zip(&b).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
        println!("N = {:3}: {}", n, sci(diff, 1));
    }
}
