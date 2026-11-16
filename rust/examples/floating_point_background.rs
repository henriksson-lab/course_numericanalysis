fn main() {
    println!("{}", 0.1 + 0.2 == 0.3);
    println!("{:.20}", 0.1 + 0.2);       // show 20 decimals of the stored value
    println!("{:e}", f64::EPSILON);      // machine epsilon; {:e} = scientific notation
    // std has no all.equal/isclose: compare with a relative tolerance yourself
    let (a, b): (f64, f64) = (0.1 + 0.2, 0.3);
    println!("{}", (a - b).abs() <= 1e-9 * a.abs().max(b.abs()));
}
