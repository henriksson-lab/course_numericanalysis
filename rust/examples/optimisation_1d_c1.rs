// C1: golden-section search, hand-written, recording the bracket width per step.
const T: [f64; 7] = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0];
const Y: [f64; 7] = [1.000, 0.804, 0.684, 0.484, 0.217, 0.094, 0.083];

fn sse(k: f64) -> f64 {
    T.iter().zip(Y.iter()).map(|(t, y)| (y - (-k * t).exp()).powi(2)).sum()
}

// returns (midpoint of final bracket, number of evaluations, bracket widths)
fn golden(f: impl Fn(f64) -> f64, mut a: f64, mut b: f64, tol: f64) -> (f64, usize, Vec<f64>) {
    let g = (5f64.sqrt() - 1.0) / 2.0;
    let (mut x1, mut x2) = (b - g * (b - a), a + g * (b - a));
    let (mut f1, mut f2) = (f(x1), f(x2));
    let mut nev = 2;
    let mut width = vec![b - a];
    while b - a > tol {
        if f1 < f2 {
            b = x2; x2 = x1; f2 = f1;
            x1 = b - g * (b - a); f1 = f(x1);
        } else {
            a = x1; x1 = x2; f1 = f2;
            x2 = a + g * (b - a); f2 = f(x2);
        }
        nev += 1;
        width.push(b - a);
    }
    ((a + b) / 2.0, nev, width)
}

fn main() {
    let (xmin, nev, width) = golden(sse, 0.0, 1.0, 1e-8);
    println!("{:.7} {}", xmin, nev);
    let pred = 2 + (1e-8f64.ln() / 0.618034f64.ln()).ceil() as usize; // predicted evaluations
    println!("{}", pred);
    let ratios: Vec<f64> = width.windows(2).map(|w| w[1] / w[0]).collect();
    let lo = ratios.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = ratios.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    println!("{:.6} {:.6}", lo, hi);
}
