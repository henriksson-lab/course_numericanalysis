// C4: the temperature step (assumed: k doubles in the warm half).

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

fn main() {
    let (dd, l, j) = (100.0_f64, 500.0_f64, 1.0_f64);
    let n = 500;
    let dx = l / n as f64;
    let r = dd / (dx * dx);
    let x: Vec<f64> = (0..=n).map(|i| i as f64 * dx).collect();
    let mut lo = vec![-r; n + 1];
    lo[0] = 0.0;
    lo[n] = -2.0 * r;
    let mut up = vec![-r; n + 1];
    up[0] = -2.0 * r;
    up[n] = 0.0;
    let mut src = vec![0.0; n + 1];
    src[0] = 2.0 * j / dx;
    // first crossing of the threshold, linear interpolation
    let x_threshold = |c: &[f64], thr: f64| -> f64 {
        let i = c.iter().position(|&v| v < thr).unwrap();
        x[i - 1] + (c[i - 1] - thr) / (c[i - 1] - c[i]) * dx
    };
    let cases: [(&str, Box<dyn Fn(f64) -> f64>); 4] = [
        ("cold", Box::new(|_| 0.01)),
        ("warm", Box::new(|_| 0.02)),
        ("posterior_warm", Box::new(move |x| if x < l / 2.0 { 0.01 } else { 0.02 })),
        ("anterior_warm", Box::new(move |x| if x < l / 2.0 { 0.02 } else { 0.01 })),
    ];
    for (nm, kf) in cases.iter() {
        let di: Vec<f64> = x.iter().map(|&xi| 2.0 * r + kf(xi)).collect();
        let c = thomas(&lo, &di, &up, &src);
        println!("{:<15} c(0) = {:.3}  c(250) = {:.4}  x(c = 0.1) = {:.1} µm", nm, c[0], c[250], x_threshold(&c, 0.1));
    }
}
