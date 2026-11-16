fn main() {
    let x = 2f64.powi(53);
    println!("{}", x + 1.0 - x);
    println!("{}", x - x + 1.0);
}
