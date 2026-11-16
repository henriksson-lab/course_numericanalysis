fn kappa_num(f: impl Fn(f64) -> f64, x: f64, h: f64) -> f64 {
    let rel_out = (f(x * (1.0 + h)) - f(x)).abs() / f(x).abs();
    rel_out / h
}

fn main() {
    let h = 1e-7;
    // |x| ... is a closure: an anonymous function, like R's function(x) ...
    println!("{}", kappa_num(|x| x - 1.0, 1.01, h));
    println!("{}", kappa_num(|x| x - 1.0, 100.0, h));
    println!("{}", kappa_num(f64::ln, 1.01, h));
    println!("{}", kappa_num(f64::ln, 100.0, h));
    println!("{}", kappa_num(f64::sqrt, 7.0, h));
    println!("{}", kappa_num(|x| 1.0 / (1.0 - x), 0.99, h));

    let k = 10.0;
    let dose = |r: f64| k * r / (1.0 - r); // the closure captures k
    for r in [0.5, 0.9, 0.99, 0.999] {
        println!("{} {} {}", r, kappa_num(dose, r, h), 1.0 / (1.0 - r));
    }
}
