// Gaussian elimination WITH partial pivoting (hand-written), then nalgebra's LU.
use nalgebra::{matrix, vector};

fn solve_pivot(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n - 1 {
        // row index p >= k of the largest |a[p][k]|
        let p = (k..n).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs())).unwrap();
        a.swap(k, p); // swap rows k and p
        b.swap(k, p);
        for i in k + 1..n {
            let m = a[i][k] / a[k][k]; // now |m| <= 1
            for j in 0..n {
                a[i][j] -= m * a[k][j];
            }
            b[i] -= m * b[k];
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let s: f64 = (i + 1..n).map(|j| a[i][j] * x[j]).sum();
        x[i] = (b[i] - s) / a[i][i];
    }
    x
}

fn rel(u: &[f64], v: &[f64]) -> f64 {
    let num: f64 = u.iter().zip(v).map(|(a, b)| (a - b).powi(2)).sum();
    let den: f64 = v.iter().map(|a| a * a).sum();
    (num / den).sqrt()
}

fn main() {
    for k in [2, 8, 14, 20] {
        let d = 10f64.powi(-k);
        let a = vec![vec![d, 1.0], vec![1.0, 1.0]];
        let x = [1.0 / (1.0 - d), (1.0 - 2.0 * d) / (1.0 - d)];
        println!("delta = 1e-{:02}  error = {:.1e}", k, rel(&solve_pivot(a, vec![1.0, 2.0]), &x));
    }

    let b3 = vec![vec![2.0, 1.0, 1.0], vec![4.0, 3.0, 3.0], vec![8.0, 7.0, 9.0]];
    println!("{:?}", solve_pivot(b3, vec![4.0, 10.0, 24.0]));

    // nalgebra: PA = LU with partial pivoting
    let bm = matrix![2.0, 1.0, 1.0;
                     4.0, 3.0, 3.0;
                     8.0, 7.0, 9.0];
    let lu = bm.lu();
    print!("L ={}U ={}", lu.l(), lu.u());
    let mut pb = bm;
    lu.p().permute_rows(&mut pb); // rows of B in pivot order
    print!("PB ={}", pb);
    println!("{}", lu.solve(&vector![4.0, 10.0, 24.0]).unwrap().transpose());
}

// Output:
// delta = 1e-02  error = 1.6e-16
// delta = 1e-08  error = 1.6e-16
// delta = 1e-14  error = 0.0e0
// delta = 1e-20  error = 0.0e0
// [1.0, 1.0, 1.0]
// L =
//   ┌                                                          ┐
//   │                  1                  0                  0 │
//   │               0.25                  1                  0 │
//   │                0.5 0.6666666666666666                  1 │
//   └                                                          ┘
//
// U =
//   ┌                                                             ┐
//   │                   8                   7                   9 │
//   │                   0               -0.75               -1.25 │
//   │                   0                   0 -0.6666666666666667 │
//   └                                                             ┘
//
// PB =
//   ┌       ┐
//   │ 8 7 9 │
//   │ 2 1 1 │
//   │ 4 3 3 │
//   └       ┘
//
//
//   ┌       ┐
//   │ 1 1 1 │
//   └       ┘
