// The circadian gene: fft, amplitude, peak time and the inverse (rustfft).
use rustfft::{num_complex::Complex, FftPlanner};

fn main() {
    let x = [7.81, 8.81, 9.05, 9.98, 9.40, 8.50, 8.15, 7.47, 6.87, 6.41, 7.15, 7.37,
             7.81, 8.09, 9.64, 9.49, 9.29, 9.03, 8.25, 7.43, 6.98, 6.73, 6.72, 6.65];
    let n = x.len();                     // 24 samples, every 2 h, T = 48 h
    let mut planner = FftPlanner::<f64>::new();
    // rustfft works in place on complex numbers: copy x into a complex buffer
    let mut xc: Vec<Complex<f64>> = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
    planner.plan_fft_forward(n).process(&mut xc);   // unnormalised: xc[0] = sum(x)
    println!("{}", xc[0].re / n as f64); // mean
    let k = 2;                           // 2 cycles per 48 h = 24-h period
    println!("{}", 2.0 * xc[k].norm() / n as f64);  // amplitude (index k, 0-based)
    let peak = (-xc[k].arg() / (2.0 * std::f64::consts::PI) * 24.0).rem_euclid(24.0);
    println!("{}", peak);                // peak time (h)
    planner.plan_fft_inverse(n).process(&mut xc);   // inverse, also unnormalised
    let err = (0..n).map(|j| (xc[j].re / n as f64 - x[j]).abs()).fold(0.0, f64::max);
    println!("{:e}", err);
}
