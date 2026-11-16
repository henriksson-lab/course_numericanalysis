fn main() {
    let mut eps: f64 = 1.0;
    while 1.0 + eps / 2.0 != 1.0 {
        eps /= 2.0;
    }
    println!("{:e}", eps);               // {:e} = scientific notation
    println!("{:e}", f64::EPSILON);

    // gap from the formula 2^(floor(log2|x|) - 52); a closure is an inline function
    let gap = |x: f64| 2f64.powf(x.abs().log2().floor() - 52.0);
    // x.next_up() is the next double above x, so next_up() - x is the true gap
    // (like np.spacing)
    for x in [1.0, 1000.0, 0.001, 3.1e9, 1e16] {
        println!("{:e}  {:e}  {:e}", x, gap(x), x.next_up() - x);
    }
    println!("{}", 1e16 + 1.0 == 1e16);
}
