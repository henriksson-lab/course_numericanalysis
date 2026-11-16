// C1: the DFT as an n x n matrix (hand-written, O(n^2)) against rustfft.
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;
use std::time::Instant;

// The DFT matrix F[j][k] = w^(jk), w = exp(-2 pi i / n)
fn dft_matrix(n: usize) -> Vec<Vec<Complex<f64>>> {
    (0..n).map(|j| (0..n).map(|k| {
        Complex::from_polar(1.0, -2.0 * PI * ((j * k) % n) as f64 / n as f64)
    }).collect()).collect()
}

fn matvec(f: &[Vec<Complex<f64>>], x: &[f64]) -> Vec<Complex<f64>> {
    f.iter().map(|row| row.iter().zip(x).map(|(a, &b)| a * b).sum()).collect()
}

fn main() {
    let f4 = dft_matrix(4);
    let x4 = matvec(&f4, &[3.0, 1.0, 1.0, 3.0]);              // check B1
    let s: Vec<String> = x4.iter()
        .map(|z| format!("{}{:+}i", z.re.round() + 0.0, z.im.round() + 0.0)).collect();
    println!("{}", s.join(" "));
    for j in 0..4 {                                            // conj(F)^T F = 4 I
        let row: Vec<f64> = (0..4).map(|k| {
            let s: Complex<f64> = (0..4).map(|m| f4[m][j].conj() * f4[m][k]).sum();
            s.re.round() + 0.0
        }).collect();
        println!("{:?}", row);
    }
    let mut planner = FftPlanner::<f64>::new();
    for n in [256, 512, 1024, 2048] {
        let x: Vec<f64> = (0..n).map(|j| (0.1 * (j * j) as f64).sin()).collect();
        let t0 = Instant::now();
        let x1 = matvec(&dft_matrix(n), &x);
        let t_mat = t0.elapsed().as_secs_f64();
        let fft = planner.plan_fft_forward(n);
        let t0 = Instant::now();
        let mut x2 = Vec::new();
        for _ in 0..100 {
            x2 = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
            fft.process(&mut x2);
        }
        let t_fft = t0.elapsed().as_secs_f64() / 100.0;
        let diff = x1.iter().zip(&x2).map(|(a, b)| (a - b).norm()).fold(0.0, f64::max);
        println!("n = {:4}  max diff {:.1e}  matrix {:.3} s  fft {:.6} s", n, diff, t_mat, t_fft);
    }
}
