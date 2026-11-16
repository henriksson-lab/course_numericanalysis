// C5: 200 right-hand sides, one at a time: refactorise, factor once, or invert.
use nalgebra::DMatrix;
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, StandardNormal};
use std::time::Instant;

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let (n, m) = (1000, 200);
    let a = DMatrix::<f64>::from_fn(n, n, |_, _| StandardNormal.sample(&mut rng));
    let bm = DMatrix::<f64>::from_fn(n, m, |_, _| StandardNormal.sample(&mut rng));

    let t0 = Instant::now();
    for j in 0..m {
        let _x = a.clone().lu().solve(&bm.column(j)).unwrap();   // factorise every time
    }
    let t_each = t0.elapsed().as_secs_f64();

    let t0 = Instant::now();
    let f = a.clone().lu();                                      // factorise once
    for j in 0..m {
        let _x = f.solve(&bm.column(j)).unwrap();
    }
    let t_once = t0.elapsed().as_secs_f64();

    let t0 = Instant::now();
    let ai = a.clone().lu().try_inverse().unwrap();              // inverse once
    for j in 0..m {
        let _x = &ai * bm.column(j);
    }
    let t_inv = t0.elapsed().as_secs_f64();
    println!("each {t_each:.3}  once {t_once:.3}  inverse {t_inv:.3}");
}
