// rust/examples/gradient_methods_c8.rs  —  cargo run --example gradient_methods_c8
// C8: multistart: argmin's Nelder-Mead and BFGS (analytic gradient) from 15 starts.
use argmin::core::{Executor, State};
use argmin::solver::linesearch::MoreThuenteLineSearch;
use argmin::solver::neldermead::NelderMead;
use argmin::solver::quasinewton::BFGS;
include!("gradient_methods_shared/defs.rs");

fn main() {
    for hill in [0.5, 1.0, 3.0] {
        for ec50 in [0.01, 0.1, 1.0, 10.0, 100.0] {
            let s: Vec<f64> = vec![0.0, 100.0, ec50, hill];
            let mut simplex = vec![s.clone()];
            for j in 0..4 { let mut v = s.clone(); v[j] += 0.1 * s[j].abs().max(1.0); simplex.push(v); }
            let nm = Executor::new(Fit, NelderMead::new(simplex).with_sd_tolerance(1e-8).unwrap())
                .configure(|st| st.max_iters(500)).run().unwrap();
            let (nmp, nmc) = (nm.state().get_best_param().unwrap().clone(), nm.state().get_best_cost());
            let eye: Vec<Vec<f64>> = (0..4).map(|i| (0..4).map(|j| (i == j) as i32 as f64).collect()).collect();
            let bf = Executor::new(Fit, BFGS::new(MoreThuenteLineSearch::new()))
                .configure(|st| st.param(s.clone()).inv_hessian(eye).max_iters(200)).run();
            let bfc = match bf { Ok(r) => format!("{:7.2}", r.state().get_best_cost()), Err(_) => "  error".into() };
            println!("EC50 {:6} hill {:3} | NM SSE {:8.2} hill {:6.2} | BFGS SSE {}", ec50, hill, nmc, nmp[3], bfc);
        }
    }
}
