// C5: hand-written Lanczos, with and without full reorthogonalisation.
use nalgebra::{DMatrix, DVector};

struct Lanczos { alpha: Vec<f64>, beta: Vec<f64>, q: DMatrix<f64> }

fn lanczos(afun: &dyn Fn(&DVector<f64>) -> DVector<f64>, n: usize, m: usize, reorth: bool) -> Lanczos {
    let mut qm = DMatrix::<f64>::zeros(n, m);
    let (mut alpha, mut beta) = (vec![0.0; m], vec![0.0; m]);
    let mut q = DVector::from_element(n, 1.0 / (n as f64).sqrt());
    let mut q_old = DVector::<f64>::zeros(n);
    let mut b = 0.0;
    for j in 0..m {
        qm.set_column(j, &q);
        let mut w = afun(&q) - &q_old * b;
        alpha[j] = w.dot(&q);
        w -= &q * alpha[j];
        if reorth {                                   // full reorthogonalisation
            let qj = qm.columns(0, j + 1);
            w -= &qj * (qj.transpose() * &w);
        }
        b = w.norm();
        beta[j] = b;
        q_old = q;
        q = w / b;
    }
    Lanczos { alpha, beta, q: qm }
}

fn ritz(l: &Lanczos, m: usize) -> Vec<f64> {          // eigenvalues of the m x m tridiagonal T
    let t = DMatrix::from_fn(m, m, |i, j| {
        if i == j { l.alpha[i] } else if i == j + 1 || j == i + 1 { l.beta[i.min(j)] } else { 0.0 }
    });
    let mut v: Vec<f64> = t.symmetric_eigen().eigenvalues.iter().cloned().collect();
    v.sort_by(|x, y| y.partial_cmp(x).unwrap());
    v
}

fn main() {
    let mut spec = vec![10.0, 9.0, 8.0];                       // 3 strong + 197 bulk
    spec.extend((0..197).map(|i| 5.0 * i as f64 / 196.0));
    let s = DVector::from_vec(spec);
    let afun = |v: &DVector<f64>| s.component_mul(v);           // A = diag(spec), never stored
    let lr = lanczos(&afun, 200, 60, true);
    let lp = lanczos(&afun, 200, 60, false);
    for m in [5, 10, 20, 60] {
        let r: Vec<String> = ritz(&lr, m)[..3].iter().map(|x| format!("{:.6}", x)).collect();
        println!("m = {}  top Ritz values: {}", m, r.join(" "));
    }
    println!("{:.6?}", &ritz(&lp, 60)[..8]);
    println!("{:.6?}", &ritz(&lr, 60)[..8]);
    let count = |l: &Lanczos| ritz(l, 60).iter().filter(|&&x| x > 7.5).count();
    println!("{} {}", count(&lp), count(&lr));
    let orth = |l: &Lanczos| (l.q.transpose() * &l.q - DMatrix::<f64>::identity(60, 60)).amax();
    println!("{:e} {:e}", orth(&lp), orth(&lr));
}
