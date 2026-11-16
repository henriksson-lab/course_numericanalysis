// C4: time the factorisation and one solve with the factors (run with --release).
use nalgebra::{DMatrix, DVector};
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, StandardNormal};
use std::time::Instant;

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    for n in [250, 500, 1000, 2000] {
        let a = DMatrix::<f64>::from_fn(n, n, |_, _| StandardNormal.sample(&mut rng));
        let b = DVector::<f64>::from_fn(n, |_, _| StandardNormal.sample(&mut rng));
        let t0 = Instant::now();
        let mut f = a.clone().lu();
        for _ in 0..4 { f = a.clone().lu(); }
        let t_fac = t0.elapsed().as_secs_f64() / 5.0;
        let t0 = Instant::now();
        let mut x = f.solve(&b).unwrap();
        for _ in 0..99 { x = f.solve(&b).unwrap(); }
        let t_sol = t0.elapsed().as_secs_f64() / 100.0;
        assert!(x[0].is_finite());
        println!("{n:5}  factor {t_fac:.4} s  solve {t_sol:.6} s  ratio {:4.0}", t_fac / t_sol);
    }
}
