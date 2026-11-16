// Gaussian elimination WITHOUT pivoting, hand-written, on A = [delta 1; 1 1].
// Matrices are Vec<Vec<f64>> (a list of rows) to keep the indexing like R.

fn solve_nopivot(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
    let n = b.len();
    for k in 0..n - 1 {
        for i in k + 1..n {
            let m = a[i][k] / a[k][k]; // multiplier
            for j in 0..n {
                a[i][j] -= m * a[k][j];
            }
            b[i] -= m * b[k];
        }
    }
    let mut x = vec![0.0; n]; // back substitution
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
    for k in [2, 5, 8, 11, 14, 17, 20] {
        let d = 10f64.powi(-k);
        let a = vec![vec![d, 1.0], vec![1.0, 1.0]];
        let x = [1.0 / (1.0 - d), (1.0 - 2.0 * d) / (1.0 - d)]; // exact solution
        let err = rel(&solve_nopivot(a, vec![1.0, 2.0]), &x);
        println!("delta = 1e-{:02}  error = {:.1e}  eps/delta = {:.1e}",
                 k, err, f64::EPSILON / d);
    }
}

// Output:
// delta = 1e-02  error = 3.3e-15  eps/delta = 2.2e-14
// delta = 1e-05  error = 2.7e-12  eps/delta = 2.2e-11
// delta = 1e-08  error = 3.5e-9  eps/delta = 2.2e-8
// delta = 1e-11  error = 5.8e-8  eps/delta = 2.2e-5
// delta = 1e-14  error = 5.7e-4  eps/delta = 2.2e-2
// delta = 1e-17  error = 7.1e-1  eps/delta = 2.2e1
// delta = 1e-20  error = 7.1e-1  eps/delta = 2.2e4
