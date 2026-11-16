// C3: Richardson table for the central difference (hand-written).
const K: f64 = 1.2;
const N0: f64 = 0.05;
const R: f64 = 0.8;

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

fn richardson(f: impl Fn(f64) -> f64, x: f64, h0: f64, levels: usize) -> Vec<Vec<f64>> {
    let mut a = vec![vec![f64::NAN; levels]; levels];
    for i in 0..levels {
        let h = h0 / 2f64.powi(i as i32);
        a[i][0] = (f(x + h) - f(x - h)) / (2.0 * h); // central difference
    }
    for j in 1..levels {
        let p = 4f64.powi(j as i32); // remove the h^(2j) term
        for i in j..levels {
            a[i][j] = (p * a[i][j - 1] - a[i - 1][j - 1]) / (p - 1.0);
        }
    }
    a
}

fn main() {
    let exact = R * od(3.0) * (1.0 - od(3.0) / K);
    for h0 in [0.4, 1e-4] {
        for row in richardson(od, 3.0, h0, 4) {
            let cells: Vec<String> = row.iter()
                .map(|v| if v.is_nan() { "      NA".into() } else { format!("{:>9.2e}", v - exact) })
                .collect();
            println!("{}", cells.join(" "));
        }
        println!();
    }
}
