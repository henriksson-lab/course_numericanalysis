// C2: the gradient forms (IMEX Euler); local accumulation times; the IMEX
// steady state does not depend on dt.

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

/// IMEX Euler: diffusion implicit (one Thomas solve), reaction f explicit.
fn imex(c: &[f64], g: &Grid, dt: f64, nsteps: usize, f: &dyn Fn(&[f64]) -> Vec<f64>) -> Vec<f64> {
    let a: Vec<f64> = g.lo.iter().map(|v| dt * v).collect();
    let b: Vec<f64> = g.di.iter().map(|v| 1.0 + dt * v).collect();
    let cc: Vec<f64> = g.up.iter().map(|v| dt * v).collect();
    let mut c = c.to_vec();
    for _ in 0..nsteps {
        let fc = f(&c);
        let rhs: Vec<f64> = (0..c.len()).map(|i| c[i] + dt * (fc[i] + g.src[i])).collect();
        c = thomas(&a, &b, &cc, &rhs);
    }
    c
}

fn main() {
    let (dd, k, l) = (100.0, 0.01, 500.0);
    let lam: f64 = (dd / k as f64).sqrt();
    let g = grid(50, dd, l, dd / lam);
    let n = g.n;
    let di_k: Vec<f64> = g.di.iter().map(|v| v + k).collect();
    let cs = thomas(&g.lo, &di_k, &g.up, &g.src); // steady state
    let f = |c: &[f64]| c.iter().map(|v| -k * v).collect::<Vec<f64>>();
    let dt = 0.5;
    let mut c = vec![0.0; n + 1];
    let mut tau = vec![0.0; n + 1];
    for _ in 0..2000 { // to t = 1000 min
        let cnew = imex(&c, &g, dt, 1, &f);
        for i in 0..=n {
            tau[i] += dt * (1.0 - (c[i] + cnew[i]) / (2.0 * cs[i])); // trapezoid rule
        }
        c = cnew;
    }
    let ix = [0, 10, 20]; // x = 0, 100, 200 µm
    println!("numeric {:?}", ix.map(|i| tau[i]));
    println!("formula {:?}", ix.map(|i| (1.0 + g.x[i] / lam) / (2.0 * k)));
    for dt in [1.0, 10.0, 50.0] {
        let c = imex(&vec![0.0; n + 1], &g, dt, (5000.0 / dt) as usize, &f);
        let e = c.iter().zip(&cs).map(|(p, q)| (p - q).abs()).fold(0.0, f64::max);
        println!("{} {:e}", dt, e);
    }
}
