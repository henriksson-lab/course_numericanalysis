// C6: two formulas, equal in exact arithmetic, give different chaotic trajectories.
const K: f64 = 100.0;

fn ricker_a(n: f64, r: f64) -> f64 { n * (r * (1.0 - n / K)).exp() }
fn ricker_b(n: f64, r: f64) -> f64 { n * r.exp() * (-r * n / K).exp() } // the same in exact arithmetic

fn main() {
    for r in [1.5, 3.0] {
        let (mut a, mut b) = (20.0, 20.0);
        let mut d = vec![0.0; 100];
        for t in 0..100 {
            a = ricker_a(a, r);
            b = ricker_b(b, r);
            d[t] = (a - b).abs();
        }
        let shown: Vec<String> = [0, 19, 39, 59, 79, 99].iter().map(|&t| format!("{:.1e}", d[t])).collect();
        println!("r = {r} : {}", shown.join(" "));
    }
}
