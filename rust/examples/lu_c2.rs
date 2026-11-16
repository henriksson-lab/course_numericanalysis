// C2: nalgebra's LU with partial pivoting: P A = L U.
use nalgebra::{dmatrix, DMatrix};

fn show(m: &DMatrix<f64>) {
    for r in m.row_iter() {
        println!("{:?}", r.iter().collect::<Vec<_>>());
    }
}

fn main() {
    let a = dmatrix![1.0, 1.0, 1.0; -13.0, -21.0, -27.0; 6.0, 18.0, 2.0];
    let b = dmatrix![1.0, 1.0; -18.4, -22.4; 11.6, 7.6];   // two mussels, one per column
    let f = a.clone().lu();
    let (l, u) = (f.l(), f.u());
    show(&l);
    show(&u);
    let mut plu = &l * &u;              // = P A: undo the row swaps to get A back
    f.p().inv_permute_rows(&mut plu);
    println!("{:e}", (plu - &a).abs().max());
    show(&f.solve(&b).expect("singular"));   // both mussels, one factorisation
}
