// C3: the T-sensor by marching in x. Implicit (backward Euler) marching checked against
// the exact plug-flow solution; flux conservation; explicit marching and its step limit.
use std::f64::consts::PI;

const W: f64 = 100.0; // channel width, um
const UBAR: f64 = 1000.0; // mean speed, um/s
const D: f64 = 1000.0; // small molecule, um^2/s

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

// backward Euler in x: one tridiagonal solve per step of dx
fn march(c: &[f64], u: &[f64], d: f64, dy: f64, dx: f64, nsteps: usize) -> Vec<f64> {
    let m = c.len();
    let a: Vec<f64> = u.iter().map(|uj| d * dx / (uj * dy * dy)).collect();
    let mut di: Vec<f64> = a.iter().map(|aj| 1.0 + 2.0 * aj).collect();
    di[0] = 1.0 + a[0]; // no-flux walls
    di[m - 1] = 1.0 + a[m - 1];
    let off: Vec<f64> = a.iter().map(|aj| -aj).collect();
    let mut c = c.to_vec();
    for _ in 0..nsteps {
        c = thomas(&off, &di, &off, &c);
    }
    c
}

fn explicit(c: &[f64], u: &[f64], d: f64, dy: f64, dx: f64, nsteps: usize) -> Vec<f64> {
    let m = c.len();
    let mut c = c.to_vec();
    for _ in 0..nsteps {
        let old = c.clone();
        for j in 0..m {
            let left = if j == 0 { old[0] } else { old[j - 1] }; // mirror ghost cells
            let right = if j == m - 1 { old[m - 1] } else { old[j + 1] };
            c[j] = old[j] + d * dx / (u[j] * dy * dy) * (left - 2.0 * old[j] + right);
        }
    }
    c
}

fn plug_exact(y: f64, x: f64, d: f64) -> f64 {
    let mut s = 0.5;
    for n in 1..=2001 {
        let nf = n as f64;
        s += 2.0 / (nf * PI) * (nf * PI / 2.0).sin() * (nf * PI * y / W).cos()
            * (-(nf * PI).powi(2) * d * x / (UBAR * W * W)).exp();
    }
    s
}

fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn grid(m: usize) -> (f64, Vec<f64>, Vec<f64>) {
    let dy = W / m as f64;
    let y: Vec<f64> = (1..=m).map(|j| (j as f64 - 0.5) * dy).collect();
    let c0 = y.iter().map(|&yj| if yj < W / 2.0 { 1.0 } else { 0.0 }).collect();
    (dy, y, c0)
}

fn maxabs(v: &[f64]) -> f64 { v.iter().fold(0.0, |m, x| m.max(x.abs())) }

fn main() {
    // (a) plug flow against the exact solution at x = 2 mm
    let (dy, y, c0) = grid(200);
    let ex: Vec<f64> = y.iter().map(|&yj| plug_exact(yj, 2000.0, D)).collect();
    for dx in [200.0, 100.0, 50.0, 25.0, 12.5] {
        let cm = march(&c0, &vec![UBAR; 200], D, dy, dx, (2000.0 / dx) as usize);
        let err = cm.iter().zip(&ex).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max);
        println!("dx = {:5.1}  max error {}", dx, sci(err, 2));
    }
    // (b) Poiseuille: flux conservation and the explicit step limits
    let (dy, y, c0) = grid(50);
    let u: Vec<f64> = y.iter().map(|yj| 6.0 * UBAR * (yj / W) * (1.0 - yj / W)).collect();
    let c2 = march(&c0, &u, D, dy, 50.0, 40);
    let flux = |c: &[f64]| u.iter().zip(c).map(|(a, b)| a * b).sum::<f64>();
    println!("flux ratio {:.12}", flux(&c2) / flux(&c0));
    let lim: Vec<f64> = u.iter().map(|uj| uj * dy * dy / (2.0 * D)).collect();
    println!("explicit dx limit: wall row {:.4} um, centre row {:.4} um",
             lim.iter().cloned().fold(f64::INFINITY, f64::min), lim.iter().cloned().fold(0.0, f64::max));
    // (c) explicit marching, dx = 1 um and 0.1 um, to x = 100 um
    println!("explicit, dx = 1:   max|c| = {}", sci(maxabs(&explicit(&c0, &u, D, dy, 1.0, 100)), 2));
    println!("explicit, dx = 0.1: max|c| = {:.3}", maxabs(&explicit(&c0, &u, D, dy, 0.1, 1000)));
}
