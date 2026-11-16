// C3: Lie and Strang splitting with exact sub-steps (nalgebra's matrix
// exponential); error against the exact solution at t = 60 min.
use nalgebra::{DMatrix, DVector};

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

/// The grid of the lecture: N = 50 intervals of 10 µm, ghost points at both ends.
/// lo, di, up = the three diagonals of -D d2/dx2; src = influx j at x = 0.
#[allow(dead_code)]
struct Grid { n: usize, dx: f64, x: Vec<f64>, lo: Vec<f64>, di: Vec<f64>, up: Vec<f64>, src: Vec<f64> }

fn grid(n: usize, dd: f64, l: f64, j: f64) -> Grid {
    let dx = l / n as f64;
    let r = dd / (dx * dx);
    let mut lo = vec![-r; n + 1];
    lo[0] = 0.0;
    lo[n] = -2.0 * r;
    let mut up = vec![-r; n + 1];
    up[0] = -2.0 * r;
    up[n] = 0.0;
    let mut src = vec![0.0; n + 1];
    src[0] = 2.0 * j / dx;
    Grid { n, dx, x: (0..=n).map(|i| i as f64 * dx).collect(), lo, di: vec![2.0 * r; n + 1], up, src }
}

fn split_err(g: &Grid, kx: &[f64], b: &[f64], c0: &[f64], hs: &[f64]) {
    let m = g.n + 1;
    let tend = 60.0;
    // A = D d2/dx2 with ghost points (minus the positive operator)
    let a = DMatrix::from_fn(m, m, |i, j| {
        if i == j { -g.di[i] } else if i == j + 1 { -g.lo[i] } else if j == i + 1 { -g.up[i] } else { 0.0 }
    });
    let mm = &a - DMatrix::from_diagonal(&DVector::from_column_slice(kx));
    let bv = DVector::from_column_slice(b);
    let cst = -mm.clone().lu().solve(&bv).unwrap();
    let cref = &cst + (&mm * tend).exp() * (DVector::from_column_slice(c0) - &cst); // exact
    let mut aug = DMatrix::zeros(m + 1, m + 1); // diffusion + influx, done exactly
    aug.view_mut((0, 0), (m, m)).copy_from(&a);
    aug.view_mut((0, m), (m, 1)).copy_from(&bv);
    for &h in hs {
        let e = (&aug * h).exp();
        let dif = |c: &DVector<f64>| -> DVector<f64> {
            let mut ce = DVector::from_element(m + 1, 1.0);
            ce.rows_mut(0, m).copy_from(c);
            (&e * ce).rows(0, m).into_owned()
        };
        let decay = |c: &DVector<f64>, t: f64| -> DVector<f64> {
            DVector::from_fn(m, |i, _| (-kx[i] * t).exp() * c[i])
        };
        let mut lie = DVector::from_column_slice(c0);
        let mut strang = lie.clone();
        for _ in 0..(tend / h).round() as usize {
            lie = decay(&dif(&lie), h);
            strang = decay(&dif(&decay(&strang, h / 2.0)), h / 2.0);
        }
        println!("dt = {:5.2}  Lie {}  Strang {}", h, sci((&lie - &cref).amax(), 2), sci((&strang - &cref).amax(), 2));
    }
}

fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let (dd, k, l) = (100.0, 0.01, 500.0);
    let g = grid(50, dd, l, dd / (dd / k as f64).sqrt());
    let m = g.n + 1;
    let di_k: Vec<f64> = g.di.iter().map(|v| v + k).collect();
    let cs = thomas(&g.lo, &di_k, &g.up, &g.src);
    let k_step: Vec<f64> = g.x.iter().map(|&x| if x < l / 2.0 { 0.01 } else { 0.02 }).collect();
    let none = vec![0.0; m];
    split_err(&g, &k_step, &g.src, &vec![0.0; m], &[20.0, 10.0, 5.0, 2.5, 1.25]); // (a)
    split_err(&g, &vec![k; m], &none, &cs, &[20.0, 5.0]);                         // (b)
    split_err(&g, &k_step, &none, &cs, &[20.0, 5.0, 1.25]);                       // (c)
}
