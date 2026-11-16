// Q3: a tiny LCG, x_{n+1} = (5 x_n + 1) mod 8, seed 0.
fn main() {
    let mut x: u64 = 0;
    for _ in 0..9 {
        print!("{} ", x);
        x = (5 * x + 1) % 8;
    }
    println!();
}
