// C1: a 1D boundary value problem in a dead-end trench (open end c = c0,
// closed end no flux by a ghost point), linear uptake k c, solved with the
// Thomas algorithm and compared with the exact cosh solution.

// Thomas algorithm (as in Structure: symmetric, banded, sparse):
// a = sub-diagonal (a[0] unused), b = diagonal, c = super-diagonal (c[n-1] unused)
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = b.len();
    let (mut bb, mut dd) = (b.to_vec(), d.to_vec());
    for i in 1..n {
        let w = a[i] / bb[i - 1];
        bb[i] -= w * c[i - 1];
        dd[i] -= w * dd[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = dd[n - 1] / bb[n - 1];
    for i in (0..n - 1).rev() { x[i] = (dd[i] - c[i] * x[i + 1]) / bb[i]; }
    x
}

fn main() {
    let (l, d, k, c0): (f64, f64, f64, f64) = (25.0, 600.0, 6.0, 1.0);   // µm, µm^2/s, 1/s, relative
    let lam = (d / k).sqrt();                          // 10 µm
    let mut errs = Vec::new();
    let mut last = 0.0;
    for n in [10usize, 20, 40, 80, 160] {
        let h = l / n as f64;
        let s = d / (h * h);
        let mut a = vec![-s; n];
        a[n - 1] = -2.0 * s;                           // ghost point c_{N+1} = c_{N-1}
        let b = vec![2.0 * s + k; n];
        let c = vec![-s; n];
        let mut rhs = vec![0.0; n];
        rhs[0] = s * c0;                               // the open end is known
        let cn = thomas(&a, &b, &c, &rhs);
        let err = (0..n).map(|i| {
            let x = h * (i + 1) as f64;
            (cn[i] - c0 * ((l - x) / lam).cosh() / (l / lam).cosh()).abs()
        }).fold(0.0, f64::max);
        errs.push(err);
        last = cn[n - 1];
    }
    println!("{:?}", errs.iter().map(|e| format!("{:.6e}", e)).collect::<Vec<_>>());
    println!("{:?}", errs.windows(2).map(|w| format!("{:.6}", w[0] / w[1])).collect::<Vec<_>>());
    println!("{} {}", last, c0 / (l / lam).cosh());    // closed end: grid vs exact
}
