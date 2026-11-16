// rust/examples/iterative_solvers_background.rs  —  cargo run --example iterative_solvers_background
// The chamber matrix in sprs. sprs has no sparse direct solver, so the
// reference solution comes from (hand-written) CG with a tight tolerance.
include!("iterative_solvers_shared/defs.rs");

fn main() {
    let m = 63;
    let (a, b) = (lap2(m), chamber_rhs(m, 12.0));
    println!("{} {}", a.rows(), a.nnz());
    let (u, _) = pcg(&a, &b, |r| r.to_vec(), 1e-12, 10_000);
    println!("{}", 1.0 - u[(m * m - 1) / 2]);       // oxygen in the centre
}
