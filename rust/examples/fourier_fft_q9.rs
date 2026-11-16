// Quiz: fft of one cycle of a cosine sampled 4 times.
use rustfft::{num_complex::Complex, FftPlanner};

fn main() {
    let mut x: Vec<Complex<f64>> = [1.0, 0.0, -1.0, 0.0].iter().map(|&v| Complex::new(v, 0.0)).collect();
    FftPlanner::<f64>::new().plan_fft_forward(4).process(&mut x);
    let s: Vec<String> = x.iter().map(|z| format!("{}{:+}i", z.re + 0.0, z.im + 0.0)).collect();
    println!("{}", s.join(" "));
    println!("{}", 2.0 * x[1].norm() / 4.0);
}
