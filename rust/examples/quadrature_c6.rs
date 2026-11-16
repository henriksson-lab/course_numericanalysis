// C6: area and height of a Gaussian chromatography peak (sigma = 3 s) from
// samples every h seconds, worst case over 40 sampling phases.
fn peak(t: f64) -> f64 {
    1000.0 * (-(t - 300.0).powi(2) / (2.0 * 9.0)).exp()
}

fn main() {
    let area_true = 1000.0 * 3.0 * (2.0 * std::f64::consts::PI).sqrt();
    for h in [9.0, 6.0, 4.5, 3.0, 1.5] {
        let (mut worst_area, mut worst_height) = (0.0f64, 0.0f64);
        for k in 0..40 {
            let off = h * k as f64 / 40.0;
            let x: Vec<f64> = (0..).map(|i| 240.0 + off + i as f64 * h).take_while(|&t| t <= 360.0 + 1e-9).collect();
            let y: Vec<f64> = x.iter().map(|&t| peak(t)).collect();
            let area = h * (y.iter().sum::<f64>() - (y[0] + y[y.len() - 1]) / 2.0);
            let ymax = y.iter().cloned().fold(f64::MIN, f64::max);
            worst_area = worst_area.max((area / area_true - 1.0).abs());
            worst_height = worst_height.min(ymax / 1000.0 - 1.0);
        }
        println!("h = {:3.1} s ({:3.1} sigma): worst area error {:.1e}, worst height error {:3.0}%",
                 h, h / 3.0, worst_area, 100.0 * worst_height);
    }
}
