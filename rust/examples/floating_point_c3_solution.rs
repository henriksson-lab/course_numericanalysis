fn main() {
    // like seq(0, 1, by = 0.1): the i-th point is i * 0.1
    let conc: Vec<f64> = (0..=10).map(|i| i as f64 * 0.1).collect();
    // position() returns Option<usize>: Some(index) or None; indices start at 0
    println!("{:?}", conc.iter().position(|&c| c == 0.3));

    let near = |a: f64, b: f64, tol: f64| (a - b).abs() <= tol * a.abs().max(b.abs()).max(1.0);
    let i = conc.iter().position(|&c| near(c, 0.3, 1e-9));
    println!("{:?}", i);
    println!("{:.20}", conc[i.unwrap()]);
}
