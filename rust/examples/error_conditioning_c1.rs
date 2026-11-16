use std::f64::consts::E; // e = 2.71828...

fn rel_err(approx: f64, truth: f64) -> f64 {
    (approx - truth).abs() / truth.abs()
}

fn digits(approx: f64, truth: f64) -> f64 {
    -rel_err(approx, truth).log10()
}

fn main() {
    let approx = [2.718, 2.7183, 19.0 / 7.0, 2.71828183];
    // map turns each number into text ({:.6e} = scientific, 6 decimals),
    // collect gathers the texts in a Vec, join puts them on one line
    let r: Vec<String> = approx.iter().map(|&a| format!("{:.6e}", rel_err(a, E))).collect();
    let d: Vec<String> = approx.iter().map(|&a| format!("{:.6}", digits(a, E))).collect();
    println!("{}", r.join("  "));
    println!("{}", d.join("  "));
}
