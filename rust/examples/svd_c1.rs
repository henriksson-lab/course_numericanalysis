// C1: the SVD of A = [[3, 0], [4, 5]], checked; and why not eigenvalues of A^T A.
use nalgebra::DMatrix;

fn show(m: &DMatrix<f64>) {                  // print a matrix row by row
    for i in 0..m.nrows() {
        let row: Vec<String> = m.row(i).iter().map(|x| format!("{:11.8}", x)).collect();
        println!("{}", row.join(" "));
    }
}

fn main() {
    let a = DMatrix::from_row_slice(2, 2, &[3.0, 0.0, 4.0, 5.0]);
    let svd = a.clone().svd(true, true);
    let d = svd.singular_values.clone();
    let mut u = svd.u.unwrap();
    let mut v = svd.v_t.unwrap().transpose();
    for j in 0..2 {                          // sign convention: first entry of each v_j > 0
        if v[(0, j)] < 0.0 {
            v.column_mut(j).neg_mut();
            u.column_mut(j).neg_mut();
        }
    }
    println!("{:?}", d.as_slice());
    show(&u);
    show(&v);
    println!("max |U^T U - I| = {:e}", (u.transpose() * &u - DMatrix::<f64>::identity(2, 2)).abs().max());
    println!("sigma^2 = {:?}", d.map(|x| x * x).as_slice());
    let ev = (a.transpose() * &a).symmetric_eigen().eigenvalues;
    println!("eigen(A^T A) = {:?}", ev.as_slice());   // not sorted

    let b = DMatrix::from_row_slice(3, 2, &[1.0, 1.0, 1.0, 1.0 + 1e-9, 1.0, 1.0 - 1e-9]);
    println!("{:?}", b.clone().singular_values().as_slice());
    let ev = (b.transpose() * &b).symmetric_eigen().eigenvalues;
    println!("{:?}", ev.map(|x: f64| x.abs().sqrt()).as_slice());
}
