// C1: steady gradient, ghost-point vs one-sided boundary: error and mass balance.

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

const DD: f64 = 100.0;
const K: f64 = 0.01;
const L: f64 = 500.0;

fn steady(n: usize, onesided: bool) -> (f64, f64) {
    let lam = (DD / K).sqrt();
    let j = DD / lam;
    let dx = L / n as f64;
    let r = DD / (dx * dx);
    let mut lo = vec![-r; n + 1];
    lo[0] = 0.0;
    lo[n] = -2.0 * r;
    let mut up = vec![-r; n + 1];
    up[0] = -2.0 * r;
    up[n] = 0.0;
    let mut di = vec![2.0 * r + K; n + 1];
    let mut rhs = vec![0.0; n + 1];
    rhs[0] = 2.0 * j / dx;
    if onesided { // (c1 - c0)/dx = -j/D
        di[0] = 1.0;
        up[0] = -1.0;
        rhs[0] = j * dx / DD;
    }
    let c = thomas(&lo, &di, &up, &rhs);
    let err = (0..=n)
        .map(|i| {
            let x = i as f64 * dx;
            (c[i] - j * lam / DD * ((L - x) / lam).cosh() / (L / lam).sinh()).abs()
        })
        .fold(0.0, f64::max);
    let mass = K * dx * (c.iter().sum::<f64>() - (c[0] + c[n]) / 2.0);
    (err, mass)
}

// scientific notation as R and Python print it: 1.05e-01 rather than Rust's 1.05e-1
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    for n in [5, 10, 20, 40, 80, 160, 320] {
        let (ge, gm) = steady(n, false);
        let (oe, om) = steady(n, true);
        println!("dx = {:6.2}  ghost: err {} mass {:.6}   one-sided: err {} mass {:.6}",
                 L / n as f64, sci(ge, 3), gm, sci(oe, 3), om);
    }
}
