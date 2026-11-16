// C1: stability regions of Euler, RK4 and backward Euler, as text maps, and the RK4 limits.
use nalgebra::Complex;

fn r_euler(z: Complex<f64>) -> Complex<f64> { 1.0 + z }
fn r_rk4(z: Complex<f64>) -> Complex<f64> { 1.0 + z + z * z / 2.0 + z.powi(3) / 6.0 + z.powi(4) / 24.0 }
fn r_beuler(z: Complex<f64>) -> Complex<f64> { 1.0 / (1.0 - z) }

// bisection for g(a) * g(b) < 0
fn bisect(g: impl Fn(f64) -> f64, mut a: f64, mut b: f64) -> f64 {
    while b - a > 1e-12 {
        let m = 0.5 * (a + b);
        if g(a) * g(m) <= 0.0 { b = m } else { a = m }
    }
    0.5 * (a + b)
}

fn main() {
    let methods: [(&str, fn(Complex<f64>) -> Complex<f64>); 3] =
        [("euler", r_euler), ("rk4", r_rk4), ("beuler", r_beuler)];
    // '#' = stable (|R(z)| <= 1); Re z from -4 to 2, Im z from 3 to -3, steps of 0.25 / 0.5
    for (name, r) in methods {
        println!("{}", name);
        for i in 0..13 {
            let im = 3.0 - 0.5 * i as f64;
            let row: String = (0..25).map(|j| {
                let z = Complex::new(-4.0 + 0.25 * j as f64, im);
                if r(z).norm() <= 1.0 { '#' } else { '.' }
            }).collect();
            println!("  {}", row);
        }
    }
    let rk4 = |x: f64| 1.0 + x + x * x / 2.0 + x.powi(3) / 6.0 + x.powi(4) / 24.0;
    println!("{}", bisect(|x| rk4(x) - 1.0, -3.0, -2.0));                         // real axis
    println!("{}", bisect(|s| r_rk4(Complex::new(0.0, s)).norm() - 1.0, 2.0, 3.0)); // imaginary axis
}
