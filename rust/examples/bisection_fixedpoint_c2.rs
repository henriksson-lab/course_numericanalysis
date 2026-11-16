// C2: brackets from a grid, both steady states, their stability, and the largest catch.
const R: f64 = 0.5;
const K: f64 = 100.0;

// harvested Ricker: f(N) = g(N) - N for a catch h
fn f(n: f64, h: f64) -> f64 {
    n * (R * (1.0 - n / K)).exp() - h - n
}
fn dg(n: f64) -> f64 {
    (R * (1.0 - n / K)).exp() * (1.0 - R * n / K) // g'(N)
}

// bisection as in C1, returning only the root (the bracket is checked by the caller)
fn bisect(f: impl Fn(f64) -> f64, mut a: f64, mut b: f64, tol: f64) -> f64 {
    let mut fa = f(a);
    for _ in 0..200 {
        let m = (a + b) / 2.0;
        let fm = f(m);
        if fm.signum() == fa.signum() { a = m; fa = fm; } else { b = m; }
        if b - a < tol * m.abs().max(1.0) { break; }
    }
    (a + b) / 2.0
}

// all roots of f on a grid: bisect every interval with a sign change
fn roots_on_grid(f: impl Fn(f64) -> f64, grid: &[f64]) -> Vec<f64> {
    let mut out = Vec::new();
    for w in grid.windows(2) {           // windows(2): consecutive pairs (grid[i], grid[i+1])
        if f(w[0]).signum() != f(w[1]).signum() {
            out.push(bisect(&f, w[0], w[1], 1e-12));
        }
    }
    out
}

fn main() {
    let grid: Vec<f64> = (0..=20).map(|i| 5.0 * i as f64).collect();
    for w in grid.windows(2) {
        if f(w[0], 10.0).signum() != f(w[1], 10.0).signum() {
            println!("bracket [{}, {}]", w[0], w[1]);
        }
    }
    let roots = roots_on_grid(|n| f(n, 10.0), &grid);
    println!("{:?}", roots);
    println!("{:?}", roots.iter().map(|&n| dg(n)).collect::<Vec<f64>>());

    let fine: Vec<f64> = (1..=1000).map(|i| 0.1 * i as f64).collect();
    for h in [10.0, 12.0, 14.0, 14.2, 14.3, 15.0] {
        let rts = roots_on_grid(|n| f(n, h), &fine);
        let txt: Vec<String> = rts.iter().map(|x| format!("{x:.3}")).collect();
        println!("H = {h:5.2}: {} steady states {}", rts.len(), txt.join(" "));
    }
}
