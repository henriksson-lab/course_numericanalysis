fn main() {
    // f64 is Rust's double; 1.0 / 3.0 needs the ".0" (1 / 3 would be integer division = 0)
    for v in [0.1, 0.25, 0.3, 1.0 / 3.0, 0.5] {
        println!("{:.20}", v);
    }

    let mut s: f64 = 0.0;                // `mut`: s is allowed to change
    for _ in 0..10 {
        s += 0.1;
    }
    println!("{}", s == 1.0);
    println!("{:.20}", s);
    println!("{}", (s - 1.0).abs() <= 1e-9); // compare with a tolerance
}
