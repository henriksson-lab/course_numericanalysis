// C4: a 2-variable LP by checking every corner (vertex enumeration), std only.
// All constraints are rows of A x <= b, including -x1 <= 0 and -x2 <= 0.
fn vertex_lp(obj: [f64; 2], a: &[[f64; 2]], b: &[f64]) -> ([f64; 2], f64, usize) {
    let (mut best, mut best_val, mut nv) = ([f64::NAN; 2], f64::NEG_INFINITY, 0);
    for i in 0..b.len() {
        for j in (i + 1)..b.len() {
            let det = a[i][0] * a[j][1] - a[i][1] * a[j][0];
            if det.abs() < 1e-12 { continue; } // parallel lines: no corner
            // Cramer's rule for the 2x2 system of rows i and j
            let v = [(b[i] * a[j][1] - a[i][1] * b[j]) / det,
                     (a[i][0] * b[j] - b[i] * a[j][0]) / det];
            let feasible = a.iter().zip(b).all(|(r, &bi)| r[0] * v[0] + r[1] * v[1] <= bi + 1e-9);
            if feasible {
                nv += 1;
                let val = obj[0] * v[0] + obj[1] * v[1];
                if val > best_val { best = v; best_val = val; }
            }
        }
    }
    (best, best_val, nv)
}

fn choose(n: u32, k: u32) -> f64 {
    (1..=k).fold(1.0, |acc, i| acc * (n - k + i) as f64 / i as f64)
}

fn main() {
    let mut a = vec![[1.0, 1.0], [4.0, 1.0], [-1.0, 0.0], [0.0, -1.0]];
    let mut b = vec![6.0, 12.0, 0.0, 0.0];
    let (x, v, nv) = vertex_lp([6.0, 3.0], &a, &b);
    println!("x = {} {}  ATP = {}  corners: {}", x[0], x[1], v, nv);
    a.push([1.0, 0.0]); // + oxygen
    b.push(1.5);
    let (x, v, nv) = vertex_lp([6.0, 3.0], &a, &b);
    println!("x = {} {}  ATP = {}  corners: {}", x[0], x[1], v, nv);
    for (n, k) in [(4, 2), (20, 10), (100, 50), (1000, 500)] {
        print!("{:.6e}  ", choose(n, k));
    }
    println!();
}
