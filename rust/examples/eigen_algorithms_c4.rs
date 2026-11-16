// C4: sensitivity of eigenvalues. Transit-compartment chain (nonsymmetric, defective)
// against a symmetric matrix. nalgebra has no eigenvectors for nonsymmetric matrices,
// so the check here is "perturb and recompute".
use nalgebra::{DMatrix, Matrix3};
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, StandardNormal};

fn sci(x: f64, digits: usize) -> String {
    let s = format!("{:.*e}", digits, x);
    let (m, e) = s.split_once('e').unwrap();
    let e: i32 = e.parse().unwrap();
    format!("{}e{}{:02}", m, if e < 0 { '-' } else { '+' }, e.abs())
}

fn sorted_sym(a: Matrix3<f64>) -> Vec<f64> {
    let mut v: Vec<f64> = a.symmetric_eigen().eigenvalues.iter().cloned().collect();
    v.sort_by(|x, y| x.partial_cmp(y).unwrap());
    v
}

// distances |lambda + 1| of all eigenvalues of k
fn dist_from_minus1(k: &DMatrix<f64>) -> Vec<f64> {
    k.complex_eigenvalues().iter().map(|z| ((z.re + 1.0).powi(2) + z.im.powi(2)).sqrt()).collect()
}

fn main() {
    let n = 10;
    let mut chain = DMatrix::<f64>::from_diagonal_element(n, n, -1.0);
    for i in 1..n { chain[(i, i - 1)] = 1.0; }           // transit compartments
    let a = Matrix3::new(4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0);
    let lam_a = sorted_sym(a);
    for p in (-14..=-2).step_by(2) {
        let e = 10f64.powi(p);
        let mut k = chain.clone();
        k[(0, n - 1)] = e;
        let mut s = a;
        s[(0, 2)] = e; s[(2, 0)] = e;
        let moved = dist_from_minus1(&k).into_iter().fold(0.0, f64::max);
        let sym: f64 = sorted_sym(s).iter().zip(&lam_a).map(|(x, y)| (x - y).abs()).fold(0.0, f64::max);
        println!("eps = {}  chain moved {}  (eps^(1/n) = {})  symmetric moved {}",
                 sci(e, 0), sci(moved, 3), sci(e.powf(1.0 / n as f64), 3), sci(sym, 3));
    }
    // the check: perturb at 1e-12 and recompute
    let mut k = chain.clone();
    k[(0, n - 1)] = 1e-10;
    let mut rng = StdRng::seed_from_u64(1);
    let kp = &k + DMatrix::from_fn(n, n, |_, _| {
        let z: f64 = StandardNormal.sample(&mut rng);
        1e-12 * z
    });
    let d = dist_from_minus1(&kp);
    let (lo, hi) = (d.iter().cloned().fold(f64::INFINITY, f64::min), d.iter().cloned().fold(0.0, f64::max));
    println!("{:.8} {:.8}", lo, hi);                    // was exactly 0.1 for every eigenvalue
}
