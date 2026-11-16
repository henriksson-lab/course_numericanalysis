// Background: the harvested Ricker stock of the lecture.
const R: f64 = 0.5; // growth rate
const K: f64 = 100.0; // carrying capacity (thousand fish)
const H: f64 = 10.0; // annual catch (thousand fish)

fn g(n: f64) -> f64 {
    n * (R * (1.0 - n / K)).exp() - H // harvested Ricker: stock next year
}
fn f(n: f64) -> f64 {
    g(n) - n // steady states are the roots of f
}

fn main() {
    let mut n = 30.0;
    for _ in 0..5 {
        n = g(n); // five years of the model = five fixed-point steps
    }
    println!("{n}");
    println!("{:?}", [0.0, 50.0, 100.0].map(f));
}
