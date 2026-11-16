// C6: flux variability analysis and knock-outs (minilp).
#![allow(dead_code, unused_imports)] // the FBA setup below is repeated in every C5-C8 file
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

// solver round-off: print values within 1e-9 of zero as 0
fn z(x: f64) -> f64 { if x.abs() < 1e-9 { 0.0 } else { x } }

fn main() {
    let (lo, hi) = bounds(10.0, 15.0);
    let opt = fba(lo, hi, BIO, Maximize).unwrap().0;
    let (mut lo_f, mut hi_f) = (lo, hi);
    lo_f[8] = opt; // fix growth at its maximum
    hi_f[8] = opt;
    for j in 0..9 {
        let mut o = [0.0; 9];
        o[j] = 1.0;
        let mn = fba(lo_f, hi_f, o, Minimize).unwrap().0;
        let mx = fba(lo_f, hi_f, o, Maximize).unwrap().0;
        println!("{:7} {:5.1} {:5.1}", RXN[j], z(mn), z(mx));
    }
    let ko: Vec<f64> = (0..9).map(|j| {
        let mut h = hi;
        h[j] = 0.0; // knock out reaction j
        fba(lo, h, BIO, Maximize).unwrap().0
    }).collect();
    show(&ko);
    let mut h = hi;
    h[2] = 0.0; // GLY and BYP1 together
    h[3] = 0.0;
    println!("{}", z(fba(lo, h, BIO, Maximize).unwrap().0));
}
