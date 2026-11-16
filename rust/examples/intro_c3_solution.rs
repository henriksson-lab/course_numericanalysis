// C3: the recurrence x_{n+1} = (10/3) x_n - x_{n-1} amplifies rounding error.
fn main() {
    let mut x = vec![0.0; 31]; // mut: the vector may be changed later
    x[0] = 1.0;
    x[1] = 1.0 / 3.0; // Rust indexes from 0, like Python
    for n in 1..30 {
        x[n + 1] = (10.0 / 3.0) * x[n] - x[n - 1];
    }
    for n in [0, 5, 10, 15, 20, 25, 30] {
        // Rust has no %g format; {:11.3e} prints 4 significant digits in scientific notation
        println!("{:2}  {:11.3e}  {:11.3e}", n, x[n], (1.0f64 / 3.0).powi(n as i32));
    }
}
