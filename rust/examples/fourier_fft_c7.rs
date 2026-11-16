// C7: the convolution theorem: circular vs linear (zero-padded) convolution, and the
// gain of a moving-average filter.
use rustfft::{num_complex::Complex, FftPlanner};

fn fft(x: &[f64]) -> Vec<Complex<f64>> {
    let mut buf: Vec<Complex<f64>> = x.iter().map(|&v| Complex::new(v, 0.0)).collect();
    FftPlanner::<f64>::new().plan_fft_forward(x.len()).process(&mut buf);
    buf
}

// circular convolution of two real vectors of the same length, via the FFT
fn cconv(y: &[f64], h: &[f64]) -> Vec<f64> {
    let n = y.len();
    let mut prod: Vec<Complex<f64>> = fft(y).iter().zip(fft(h)).map(|(a, b)| a * b).collect();
    FftPlanner::<f64>::new().plan_fft_inverse(n).process(&mut prod);
    prod.iter().map(|z| z.re / n as f64).collect()
}

fn main() {
    // (a) B5 by computer
    let y = [1.0, 2.0, 3.0, 4.0];
    let h = [0.5, 0.5, 0.0, 0.0];
    println!("{:?}", cconv(&y, &h));                       // circular convolution
    let s: Vec<String> = fft(&h).iter().map(|z| format!("{}{:+}i", z.re, z.im + 0.0)).collect();
    println!("{}", s.join(" "));                           // gain of the 2-point average
    // (b) linear convolution by zero-padding to length 4 + 2 - 1 = 5
    println!("{:?}", cconv(&[1.0, 2.0, 3.0, 4.0, 0.0], &[0.5, 0.5, 0.0, 0.0, 0.0]));
    // direct sum: (h * y)_j = sum_m h_m y_(j-m)
    let hh = [0.5, 0.5];
    let direct: Vec<f64> = (0..5).map(|j| (0..2).filter(|&m| j >= m && j - m < 4)
        .map(|m| hh[m] * y[j - m]).sum()).collect();
    println!("{:?}", direct);
    // (c) the gain of a moving average of width m = 8 on n = 240 samples
    let (n, m) = (240, 8);
    let hb: Vec<f64> = (0..n).map(|j| if j < m { 1.0 / m as f64 } else { 0.0 }).collect();
    let gain: Vec<f64> = fft(&hb).iter().map(|z| z.norm()).collect();
    for k in [0, 5, 10, 15, 20, 30, 40, 45, 50, 60] {
        println!("k = {:2}  period {:5.1} h  gain {:.3}", k, 120.0 / k as f64, gain[k]);
    }
}
