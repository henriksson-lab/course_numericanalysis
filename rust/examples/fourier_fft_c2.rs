// C2: a hand-written recursive FFT (powers of 2), and rustfft on awkward n.
use rustfft::{num_complex::Complex, FftPlanner};
use std::f64::consts::PI;
use std::time::Instant;

fn my_fft(x: &[Complex<f64>]) -> Vec<Complex<f64>> {
    let n = x.len();
    if n == 1 {
        return x.to_vec();
    }
    let even: Vec<Complex<f64>> = x.iter().step_by(2).cloned().collect();          // 0, 2, 4, ...
    let odd: Vec<Complex<f64>> = x.iter().skip(1).step_by(2).cloned().collect();   // 1, 3, 5, ...
    let (e, o) = (my_fft(&even), my_fft(&odd));
    let mut out = vec![Complex::new(0.0, 0.0); n];
    for k in 0..n / 2 {
        let wo = Complex::from_polar(1.0, -2.0 * PI * k as f64 / n as f64) * o[k];
        out[k] = e[k] + wo;
        out[k + n / 2] = e[k] - wo;
    }
    out
}

fn signal(n: usize) -> Vec<Complex<f64>> {
    (0..n).map(|j| Complex::new((0.1 * (j as f64).powi(2)).sin(), 0.0)).collect()
}

fn main() {
    let x4: Vec<Complex<f64>> = [3.0, 1.0, 1.0, 3.0].iter().map(|&v| Complex::new(v, 0.0)).collect();
    let s: Vec<String> = my_fft(&x4).iter()
        .map(|z| format!("{}{:+}i", z.re + 0.0, z.im + 0.0)).collect();
    println!("{}", s.join(" "));                                   // check B2
    let mut planner = FftPlanner::<f64>::new();
    let x = signal(1024);
    let mut lib = x.clone();
    planner.plan_fft_forward(1024).process(&mut lib);
    let diff = my_fft(&x).iter().zip(&lib).map(|(a, b)| (a - b).norm()).fold(0.0, f64::max);
    println!("{:e}", diff);
    for n in [65536, 65537, 65610] {                               // 65537 is prime
        let mut x = signal(n);
        let t0 = Instant::now();
        planner.plan_fft_forward(n).process(&mut x);               // planning included
        println!("n = {}  fft {:.3} s", n, t0.elapsed().as_secs_f64());
    }
}
