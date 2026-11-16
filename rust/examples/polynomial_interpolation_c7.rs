// C7: a gene-expression time course with uneven sampling times.
fn bary_weights(x: &[f64]) -> Vec<f64> {
    (0..x.len())
        .map(|j| 1.0 / (0..x.len()).filter(|&k| k != j).map(|k| x[j] - x[k]).product::<f64>())
        .collect()
}
fn bary_eval(x: &[f64], y: &[f64], w: &[f64], z: f64) -> f64 {
    let (mut num, mut den) = (0.0, 0.0);
    for j in 0..x.len() {
        let d = z - x[j];
        if d == 0.0 {
            return y[j]; // z is a node
        }
        num += w[j] * y[j] / d;
        den += w[j] / d;
    }
    num / den
}
fn interp(x: &[f64], y: &[f64], zs: &[f64]) -> Vec<f64> {
    let w = bary_weights(x);
    zs.iter().map(|&z| bary_eval(x, y, &w, z)).collect()
}
fn linear(x: &[f64], y: &[f64], z: f64) -> f64 { // like approx() / np.interp
    let i = (1..x.len()).find(|&i| z <= x[i]).unwrap();
    y[i - 1] + (y[i] - y[i - 1]) * (z - x[i - 1]) / (x[i] - x[i - 1])
}

fn main() {
    let g = |t: f64| 1.0 + 6.0 * ((-0.25 * t).exp() - (-1.5 * t).exp()); // fold change
    let ts = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 24.0]; // sampling times (h)
    let y: Vec<f64> = ts.iter().map(|&t| g(t)).collect();
    let tq = [3.0, 6.0, 10.0, 18.0]; // (a)
    let pq = interp(&ts, &y, &tq);
    for i in 0..4 {
        println!("t = {:2}  poly {:9.3}  linear {:6.3}  truth {:6.3}",
                 tq[i], pq[i], linear(&ts, &y, tq[i]), g(tq[i]));
    }
    let tt: Vec<f64> = (0..=2400).map(|i| i as f64 * 0.01).collect(); // (b)
    let p = interp(&ts, &y, &tt);
    let err: Vec<f64> = tt.iter().zip(&p).map(|(&t, &v)| (v - g(t)).abs()).collect();
    let imax = (0..err.len()).max_by(|&i, &j| err[i].total_cmp(&err[j])).unwrap();
    let pmin = p.iter().cloned().fold(f64::INFINITY, f64::min);
    let pmax = p.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("range {pmin:.1} to {pmax:.1}, max error {:.1} at {:.2} h", err[imax], tt[imax]);
    for i in 1..7 { // (c) leave one out
        let (xs, ys): (Vec<f64>, Vec<f64>) = (0..8).filter(|&k| k != i).map(|k| (ts[k], y[k])).unzip();
        println!("leave out {:4} h: poly {:+9.3}  linear {:+6.3}", ts[i],
                 interp(&xs, &ys, &[ts[i]])[0] - y[i], linear(&xs, &ys, ts[i]) - y[i]);
    }
}
