// C2: error of (e^h - 1)/h as an estimate of the slope 1, for h = 1e-1 ... 1e-12.
fn main() {
    // 1..=12 is the range 1, 2, ..., 12 (the = includes the end)
    let h: Vec<f64> = (1..=12).map(|k| 10f64.powi(-k)).collect();
    let err: Vec<f64> = h.iter().map(|hi| ((hi.exp() - 1.0) / hi - 1.0).abs()).collect();
    for (hi, ei) in h.iter().zip(&err) {
        println!("{:>6.0e}  {:.3e}", hi, ei); // {:>6.0e}: right-aligned in 6 characters
    }
    // index of the smallest error (like which.min / argmin)
    let best = (0..h.len()).min_by(|&i, &j| err[i].total_cmp(&err[j])).unwrap();
    println!("best h = {:.0e}", h[best]);
}
