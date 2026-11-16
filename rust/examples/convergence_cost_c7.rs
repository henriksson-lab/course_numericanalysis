// C7: memory of dense and sparse matrices.
// Sparse storage with sprs (compressed sparse column, 4-byte u32 indices as in R/scipy).
// The random matrix differs from R and Python, but the sizes do not.
use rand::{rngs::StdRng, seq::index, SeedableRng};
use rand_distr::{Distribution, Poisson};
use sprs::CsMatI;
use std::mem::size_of;

fn main() {
    let x = vec![0.0_f64; 1000 * 1000]; // a 1000 x 1000 matrix of zeros, stored as one Vec
    println!("{}", x.len() * size_of::<f64>());
    println!("{} {}", 20000.0 * 100000.0 * 8.0 / 1e9, 20000.0_f64.powi(2) * 8.0 / 1e9); // GB

    // genes x cells, 10 % non-zeros: in each cell, 200 random genes get count Poisson(2) + 1
    let (genes, cells, per_cell) = (2000, 10000, 200);
    let mut rng = StdRng::seed_from_u64(1);
    let pois = Poisson::new(2.0).unwrap();
    let mut indptr = vec![0u32];
    let (mut indices, mut data) = (Vec::new(), Vec::new());
    for _ in 0..cells {
        let mut rows = index::sample(&mut rng, genes, per_cell).into_vec();
        rows.sort(); // CSC wants sorted row indices within a column
        for r in rows {
            indices.push(r as u32);
            data.push(pois.sample(&mut rng) + 1.0);
        }
        indptr.push(indices.len() as u32);
    }
    let s = CsMatI::new_csc((genes, cells), indptr, indices, data);
    let sparse_bytes = s.data().len() * size_of::<f64>()
        + s.indices().len() * size_of::<u32>()
        + s.indptr().raw_storage().len() * size_of::<u32>();
    let dense_bytes = genes * cells * size_of::<f64>(); // what s.to_dense() would take
    println!("{} {}", dense_bytes as f64 / 1e6, sparse_bytes as f64 / 1e6); // MB

    // centre each gene: every zero becomes -mean, so the result must be dense
    let mut row_sum = vec![0.0; genes];
    for (&v, (r, _c)) in s.iter() {
        row_sum[r as usize] += v; // u32 index -> usize to index a Vec
    }
    let mut sc = vec![0.0_f64; genes * cells]; // dense, column-major
    for c in 0..cells {
        for r in 0..genes {
            sc[c * genes + r] = -row_sum[r] / cells as f64;
        }
    }
    for (&v, (r, c)) in s.iter() {
        sc[c as usize * genes + r as usize] += v;
    }
    println!("{}", sc.len() * size_of::<f64>() / 1_000_000); // MB
}
