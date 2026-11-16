// C6: BCA protein assay. Polynomial calibration in ug/mL and mg/mL; choose the degree.
use nalgebra::{DMatrix, DVector};

// scientific notation as R and Python print it: 2.2e+01 rather than Rust's 2.2e1
fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn vander(x: &[f64], d: usize) -> DMatrix<f64> {
    DMatrix::from_fn(x.len(), d + 1, |i, j| x[i].powi(j as i32))
}

fn qr_solve(a: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    let qr = a.clone().qr();
    let mut qtb = b.clone();
    qr.q_tr_mul(&mut qtb);
    qr.r().solve_upper_triangular(&qtb.rows(0, a.ncols()).into_owned()).unwrap()
}

fn ne_solve(a: &DMatrix<f64>, b: &DVector<f64>) -> DVector<f64> {
    (a.transpose() * a).lu().solve(&(a.transpose() * b)).unwrap()
}

fn main() {
    let conc: Vec<f64> = [0., 25., 125., 250., 500., 750., 1000., 1500., 2000.]
        .iter().flat_map(|&c| [c; 2]).collect();
    let a562 = DVector::from_vec(vec![0.087, 0.103, 0.150, 0.109, 0.226, 0.229, 0.358,
        0.344, 0.601, 0.569, 0.780, 0.789, 0.951, 0.942, 1.300, 1.238, 1.545, 1.532]);
    let mg: Vec<f64> = conc.iter().map(|c| c / 1000.0).collect();

    let v = vander(&conc, 3);                                  // cubic, ug/mL
    let show = |x: DVector<f64>| x.iter().map(|&c| sci(c, 6)).collect::<Vec<_>>().join(" ");
    println!("{}", show(ne_solve(&v, &a562)));                 // nalgebra does not check: no error
    println!("{}", show(qr_solve(&v, &a562)));
    let vs = vander(&mg, 3);                                   // cubic, mg/mL
    let show8 = |x: DVector<f64>| x.iter().map(|c| format!("{:.8}", c)).collect::<Vec<_>>().join(" ");
    println!("{}", show8(ne_solve(&vs, &a562)));
    println!("{}", show8(qr_solve(&vs, &a562)));

    for d in 1..=3 {                                           // which degree? residuals
        let a = vander(&mg, d);
        let r = &a562 - &a * qr_solve(&a, &a562);
        let signs: String = r.iter().map(|&ri| if ri > 0.0 { '+' } else { '-' }).collect();
        println!("{} {:.4} {}", d, r.norm_squared(), signs);
    }
    let c = qr_solve(&vander(&mg, 2), &a562);                  // quadratic, mg/mL
    let y = 0.70;                                              // unknown sample
    let s = (-c[1] + (c[1] * c[1] - 4.0 * c[2] * (c[0] - y)).sqrt()) / (2.0 * c[2]);
    println!("{}", s * 1000.0);                                // ug/mL
}
