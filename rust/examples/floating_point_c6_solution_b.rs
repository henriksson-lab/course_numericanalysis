fn logsumexp(a: &[f64]) -> f64 {
    let m = a.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    m + a.iter().map(|x| (x - m).exp()).sum::<f64>().ln()
}

fn main() {
    let (n, k, err) = (2000.0, 1010.0, 0.01);       // reads, alt reads, error rate
    let q: [f64; 3] = [err, 0.5, 1.0 - err];                 // P(alt read | genotype): AA, AB, BB
    let prior: [f64; 3] = [0.25, 0.5, 0.25];
    let a: Vec<f64> = (0..3)
        .map(|g| k * q[g].ln() + (n - k) * (1.0 - q[g]).ln() + prior[g].ln())
        .collect();
    println!("{:?}", a);
    println!("{}", a.iter().map(|x| x.exp()).sum::<f64>().ln());
    let lse = logsumexp(&a);
    println!("{}", lse);
    let post: Vec<f64> = a.iter().map(|x| (x - lse).exp()).collect(); // posterior
    println!("{:?}", post);
}
