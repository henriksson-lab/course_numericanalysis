// Q5: which reaction does the discrete inverse transform pick?
fn main() {
    let a = [2.0, 0.0, 3.0]; // propensities of reactions 1, 2, 3
    let a0: f64 = a.iter().sum();
    let cum: Vec<f64> = a.iter().scan(0.0, |s, &x| { *s += x; Some(*s) }).collect();
    for u in [0.3, 0.4, 0.5, 0.9] {
        let j = cum.partition_point(|&c| c <= u * a0) + 1; // like findInterval + 1
        println!("{} -> {}", u, j);
    }
}
