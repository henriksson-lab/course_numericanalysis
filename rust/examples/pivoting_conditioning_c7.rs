// Detector noise in unmixed cells: 1000 cells, 1 % noise, overlap 0.5 vs 0.95.
// Uses nalgebra, rand and rand_distr; the random numbers differ from R and Python.
use nalgebra::{matrix, DMatrix, Vector2};
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, Normal};

fn main() {
    let mut rng = StdRng::seed_from_u64(2);
    let normal = Normal::new(0.0, 0.01).unwrap();
    let noise = DMatrix::from_fn(2, 1000, |_, _| normal.sample(&mut rng)); // 2 detectors x 1000 cells
    for c in [0.5, 0.95] {
        let a = matrix![1.0, c;
                        c,   1.0];
        let clean = a * Vector2::new(1.0, 1.0);
        let mut b = noise.clone();
        for mut col in b.column_iter_mut() {
            col += clean; // readings, every cell
        }
        let a_dyn = DMatrix::from_row_slice(2, 2, &[1.0, c, c, 1.0]); // same A, any size type
        let x = a_dyn.clone().lu().solve(&b).unwrap(); // unmix all cells at once
        let d1 = x.row(0);
        let mean = d1.mean();
        let sd = (d1.iter().map(|v: &f64| (v - mean).powi(2)).sum::<f64>() / 999.0).sqrt();
        let inv = a.try_inverse().unwrap();
        let pred = 0.01 * inv.row(0).norm();
        let s = a.singular_values();
        println!("c = {:.2}  kappa = {:4.1}  sd(dye 1) = {:.4}  predicted = {:.4}",
                 c, s.max() / s.min(), sd, pred);
    }
}

// Output:
// c = 0.50  kappa =  3.0  sd(dye 1) = 0.0158  predicted = 0.0149
// c = 0.95  kappa = 39.0  sd(dye 1) = 0.1479  predicted = 0.1415
