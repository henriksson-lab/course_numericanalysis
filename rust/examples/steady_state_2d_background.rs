// Background: the 2D five-point matrix for a 3 x 2 grid, built from triplets.
// Unknown k = i + nx * j (0-based), so x (index i) runs fastest, as in R.
use sprs::{CsMat, TriMat};

fn laplacian_2d(nx: usize, ny: usize) -> CsMat<f64> {
    let idx = |i: usize, j: usize| i + nx * j;
    let mut t = TriMat::new((nx * ny, nx * ny));
    for j in 0..ny {
        for i in 0..nx {
            let k = idx(i, j);
            t.add_triplet(k, k, 4.0);
            if i > 0 { t.add_triplet(k, idx(i - 1, j), -1.0); }      // x-neighbours
            if i + 1 < nx { t.add_triplet(k, idx(i + 1, j), -1.0); }
            if j > 0 { t.add_triplet(k, idx(i, j - 1), -1.0); }      // y-neighbours
            if j + 1 < ny { t.add_triplet(k, idx(i, j + 1), -1.0); }
        }
    }
    t.to_csr()
}

fn main() {
    let (nx, ny) = (3, 2);
    let a = laplacian_2d(nx, ny);
    let dense = a.to_dense();
    for r in dense.rows() {
        println!("{:?}", r.iter().map(|v| *v as i32).collect::<Vec<_>>());
    }
    // a grid value c[i][j] at (x_i, y_j) goes to position i + nx * j
    let c = [[1.0, 4.0], [2.0, 5.0], [3.0, 6.0]];          // c[i][j]
    let mut v = vec![0.0; nx * ny];
    for j in 0..ny { for i in 0..nx { v[i + nx * j] = c[i][j]; } }
    println!("{:?}", v);
}
