// C7: the same concentration twice (ug/mL and mg/mL): a zero on the diagonal of R.
use nalgebra::DMatrix;

// scientific notation as R and Python print it: 2.2e+01 rather than Rust's 2.2e1
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn main() {
    let conc: Vec<f64> = [0., 25., 125., 250., 500., 750., 1000., 1500., 2000.]
        .iter().flat_map(|&c| [c; 2]).collect();
    let x = DMatrix::from_fn(conc.len(), 3, |i, j| match j {
        0 => 1.0,
        1 => conc[i],
        _ => conc[i] / 1000.0,
    });
    let r = x.clone().qr().r();
    let diag: Vec<String> = r.diagonal().iter().map(|&d| sci(d, 6)).collect();
    println!("{}", diag.join(" "));
    // numerical rank: diagonal entries above a relative tolerance (as R's tol = 1e-7)
    let rmax = r.diagonal().amax();
    let rank = r.diagonal().iter().filter(|d| d.abs() > 1e-7 * rmax).count();
    println!("{}", rank);
    println!("{}", x.rank(1e-7 * x.singular_values().max()));    // via the SVD
}
