// C3: linear (bisection) vs quadratic (Newton) convergence to sqrt(2).
fn main() {
    let nstep = 60;
    let (mut lo, mut hi, mut x) = (1.0_f64, 2.0_f64, 1.0_f64); // `mut`: these change
    let mut e_bis = vec![0.0; nstep];
    let mut e_new = vec![0.0; nstep];
    for k in 0..nstep {
        let mid = (lo + hi) / 2.0;
        if mid * mid > 2.0 { hi = mid } else { lo = mid }
        e_bis[k] = hi - lo; // bracket width: a guaranteed error bound
        x = (x + 2.0 / x) / 2.0;
        e_new[k] = (x - 2f64.sqrt()).abs();
    }
    // first step below 1e-12 (counting from 1)
    println!("{}", e_bis.iter().position(|&e| e < 1e-12).unwrap() + 1);
    println!("{}", e_new.iter().position(|&e| e < 1e-12).unwrap() + 1);
    let lin: Vec<f64> = (0..5).map(|k| e_bis[k + 1] / e_bis[k]).collect(); // linear: e[k+1] / e[k]
    let quad: Vec<f64> = (0..3).map(|k| e_new[k + 1] / e_new[k].powi(2)).collect(); // quadratic: e[k+1] / e[k]^2
    let digits: Vec<f64> = e_new[0..5].iter().map(|e| -e.log10()).collect(); // correct digits, Newton
    println!("{lin:?}");
    println!("{quad:?}");
    println!("{digits:?}");
}
