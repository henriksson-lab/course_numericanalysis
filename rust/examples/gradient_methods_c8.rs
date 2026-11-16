// rust/examples/gradient_methods_c8.rs  —  cargo run --example gradient_methods_c8
// C8: multistart: argmin's Nelder-Mead, and the BFGS of C5 with the analytic
// gradient, from 15 starts. argmin's Nelder-Mead stops with an internal error
// when the cost is NaN (EC50 < 0), so this cost returns +Inf there instead.
// (argmin's BFGS with its More-Thuente line search also stops with a NaN
// error from most of these starts, as in C1(c).)
use argmin::core::{Executor, State};
use argmin::solver::neldermead::NelderMead;
include!("gradient_methods_shared/defs.rs");

struct SafeFit;
impl CostFunction for SafeFit {
    type Param = Vec<f64>;
    type Output = f64;
    fn cost(&self, p: &Vec<f64>) -> Result<f64, Error> {
        let v = sse(p);
        Ok(if v.is_nan() { f64::INFINITY } else { v })
    }
}

fn bfgs(mut x: Vec<f64>) -> f64 {           // as in C5, at most 60 halvings
    let n = x.len();
    let mut b = vec![vec![0.0; n]; n];
    for i in 0..n { b[i][i] = 1.0; }
    let mut gk = grad_sse(&x);
    for _ in 0..200 {
        if norm(&gk) < 1e-4 { break; }
        let s: Vec<f64> = (0..n).map(|i| -dot(&b[i], &gk)).collect();
        let (fx, gs) = (sse(&x), dot(&gk, &s));
        let (mut a, mut halvings) = (1.0, 0);
        while !(sse(&axpy(&x, a, &s)) <= fx + 1e-4 * a * gs) && halvings < 60 { a /= 2.0; halvings += 1; }
        if halvings == 60 { break; }
        let s: Vec<f64> = s.iter().map(|v| a * v).collect();
        x = axpy(&x, 1.0, &s);
        let gn = grad_sse(&x);
        let yv: Vec<f64> = gn.iter().zip(&gk).map(|(a, b)| a - b).collect();
        gk = gn;
        let sy = dot(&yv, &s);
        if sy > 0.0 {
            let v: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j|
                (i == j) as i32 as f64 - yv[i] * s[j] / sy).collect()).collect();
            let bv: Vec<Vec<f64>> = (0..n).map(|i| (0..n).map(|j|
                (0..n).map(|l| b[i][l] * v[l][j]).sum()).collect()).collect();
            b = (0..n).map(|i| (0..n).map(|j|
                (0..n).map(|l| v[l][i] * bv[l][j]).sum::<f64>() + s[i] * s[j] / sy).collect()).collect();
        }
    }
    sse(&x)
}

fn main() {
    for hill in [0.5, 1.0, 3.0f64] {
        for ec50 in [0.01, 0.1, 1.0, 10.0, 100.0f64] {
            let s: Vec<f64> = vec![0.0, 100.0, ec50, hill];
            let mut simplex = vec![s.clone()];
            for j in 0..4 { let mut v = s.clone(); v[j] += 0.1 * s[j].abs().max(1.0); simplex.push(v); }
            let nm = Executor::new(SafeFit, NelderMead::new(simplex).with_sd_tolerance(1e-8).unwrap())
                .configure(|st| st.max_iters(500)).run().unwrap();
            let st = nm.state();
            println!("EC50 {:6} hill {:3} | NM SSE {:8.2} hill {:6.2} iter {:3} | BFGS SSE {:7.2}",
                     ec50, hill, st.get_best_cost(), st.get_best_param().unwrap()[3], st.get_iter(), bfgs(s));
        }
    }
}
