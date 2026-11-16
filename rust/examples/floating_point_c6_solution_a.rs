// &[f64] is a slice: a view of a list of numbers
fn logsumexp(a: &[f64]) -> f64 {
    let m = a.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;                   // all terms impossible
    }
    m + a.iter().map(|x| (x - m).exp()).sum::<f64>().ln()
}

fn main() {
    println!("{}", logsumexp(&[-500.0, -502.0]));
    println!("{}", logsumexp(&[2.0, 2.0, 2.0]));
    println!("{}", logsumexp(&[-1200.0, -1198.0, -1210.0]));
    println!("{}", logsumexp(&[f64::NEG_INFINITY, f64::NEG_INFINITY]));
}
