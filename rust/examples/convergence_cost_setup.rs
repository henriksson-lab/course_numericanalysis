// Setup for C1 and C2: the OD curve of the lecture.
// (The C1 and C2 files repeat these lines, so each runs on its own.)
const K: f64 = 1.2; // OD units
const N0: f64 = 0.05; // OD units
const R: f64 = 0.8; // per hour

fn od(t: f64) -> f64 {
    K / (1.0 + (K / N0 - 1.0) * (-R * t).exp())
}

fn main() {
    let t0 = 3.0;
    let exact = R * od(t0) * (1.0 - od(t0) / K); // logistic: N' = r N (1 - N/K)
    println!("{exact}"); // check: the exact derivative at t0
}
