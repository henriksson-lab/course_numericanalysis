fn main() {
    let x: f64 = 1e8;
    println!("{}", (x * x + 1.0).sqrt() - x);
    println!("{:e}", 1.0 / ((x * x + 1.0).sqrt() + x)); // {:e}: scientific notation
}
