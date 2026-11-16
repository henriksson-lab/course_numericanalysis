// C5: flux balance analysis of the toy network with minilp; rank of S with nalgebra.
#![allow(dead_code, unused_imports)] // the FBA setup below is repeated in every C5-C8 file
use nalgebra::DMatrix;
use minilp::{ComparisonOp, Error, OptimizationDirection, Problem};

const RXN: [&str; 9] = ["EX_glc", "EX_o2", "GLY", "BYP1", "BYP2", "RES", "FER", "EX_ac", "BIO"];
const S: [[f64; 9]; 6] = [
    [1., 0., -1., -1., 0., 0., 0., 0., 0.],   // G
    [0., 1., 0., 0., 0., -2., 0., 0., 0.],    // O2
    [0., 0., 0., 1., -1., 0., 0., 0., 0.],    // X
    [0., 0., 2., 0., 2., -1., -1., 0., -1.],  // P
    [0., 0., 2., 0., 2., 11., 1., 0., -5.],   // ATP
    [0., 0., 0., 0., 0., 0., 1., -1., 0.],    // Ac
];

// Optimise obj . v subject to S v = 0 and lo <= v <= hi (knock-out: hi = 0).
fn fba(lo: [f64; 9], hi: [f64; 9], obj: [f64; 9], dir: OptimizationDirection)
       -> Result<(f64, Vec<f64>), Error> {
    let mut p = Problem::new(dir);
    let v: Vec<_> = (0..9).map(|j| p.add_var(obj[j], (lo[j], hi[j]))).collect();
    for row in S.iter() {
        let terms: Vec<_> = (0..9).filter(|&j| row[j] != 0.0).map(|j| (v[j], row[j])).collect();
        p.add_constraint(terms, ComparisonOp::Eq, 0.0);
    }
    let sol = p.solve()?;
    Ok((sol.objective(), v.iter().map(|&x| sol[x]).collect()))
}

// default bounds: 0 <= v, glucose <= glc, oxygen <= o2, the rest unbounded above
fn bounds(glc: f64, o2: f64) -> ([f64; 9], [f64; 9]) {
    let mut hi = [f64::INFINITY; 9];
    hi[0] = glc;
    hi[1] = o2;
    ([0.0; 9], hi)
}

const BIO: [f64; 9] = [0., 0., 0., 0., 0., 0., 0., 0., 1.]; // objective: growth
use OptimizationDirection::{Maximize, Minimize};

fn show(v: &[f64]) {
    println!("{}", RXN.iter().map(|n| format!("{:>7}", n)).collect::<Vec<_>>().join(" "));
    println!("{}", v.iter().map(|x| format!("{:7.3}", x + 0.0)).collect::<Vec<_>>().join(" "));
}

fn main() {
    let sm = DMatrix::from_fn(6, 9, |i, j| S[i][j]);
    let r = sm.rank(1e-10);
    println!("{}", r); // independent mass balances
    println!("{}", 9 - r); // dimension of the null space of S
    let (lo, hi) = bounds(10.0, 15.0);
    let (growth, v) = fba(lo, hi, BIO, Maximize).unwrap();
    println!("{}", growth);
    show(&v);
    let resid = (0..6).map(|i| (0..9).map(|j| S[i][j] * v[j]).sum::<f64>().abs())
        .fold(0.0, f64::max);
    println!("{:e}", resid); // steady state check
}
