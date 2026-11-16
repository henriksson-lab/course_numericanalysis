// C2: check the interpolation error formula for the degree-4 interpolant.
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

// polynomials in u as coefficient vectors [a0, a1, a2, ...]
fn pmul(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    for i in 0..a.len() {
        for j in 0..b.len() {
            out[i + j] += a[i] * b[j];
        }
    }
    out
}
fn pder(a: &[f64]) -> Vec<f64> {
    (1..a.len()).map(|i| a[i] * i as f64).collect()
}
fn peval(a: &[f64], x: f64) -> f64 {
    a.iter().rev().fold(0.0, |s, &c| s * x + c) // Horner
}

fn main() {
    let g = [0.0, 1.0, -1.0]; // u - u^2
    let mut p_list = vec![g.to_vec()]; // P1
    for k in 1..5 {
        let next = pmul(&pder(&p_list[k - 1]), &g);
        p_list.push(next);
    }
    let deriv_od = |t: f64, k: usize| K * R.powi(k as i32) * peval(&p_list[k - 1], od(t) / K);
    println!("{}", deriv_od(3.0, 2)); // f''(3)
    let tt: Vec<f64> = (0..=8000).map(|i| i as f64 * 0.001).collect();
    let m2 = tt.iter().map(|&t| deriv_od(t, 2).abs()).fold(0.0, f64::max);
    let m5 = tt.iter().map(|&t| deriv_od(t, 5).abs()).fold(0.0, f64::max);
    println!("{m2} {m5}");
    // degree-4 interpolant through 0, 2, 4, 6, 8 (Newton form from C1)
    let t = [0.0, 2.0, 4.0, 6.0, 8.0];
    let mut a: Vec<f64> = t.iter().map(|&ti| od(ti)).collect();
    for k in 1..5 {
        for i in (k..5).rev() {
            a[i] = (a[i] - a[i - 1]) / (t[i] - t[i - k]);
        }
    }
    let p = |x: f64| (0..4).rev().fold(a[4], |s, k| a[k] + (x - t[k]) * s);
    let w = |x: f64| t.iter().map(|&ti| x - ti).product::<f64>();
    for tq in [1.0, 3.0] {
        println!("t = {tq}: actual {:.4}  bound {:.4}", (od(tq) - p(tq)).abs(),
                 m5 / 120.0 * w(tq).abs());
    }
    let actual: Vec<f64> = tt.iter().map(|&x| (od(x) - p(x)).abs()).collect();
    let bound: Vec<f64> = tt.iter().map(|&x| m5 / 120.0 * w(x).abs()).collect();
    println!("{}", actual.iter().zip(&bound).all(|(a, b)| a <= b));
    let ratio = actual.iter().zip(&bound).filter(|(_, b)| **b > 0.0)
        .map(|(a, b)| a / b).fold(0.0, f64::max);
    println!("{ratio}");
}
