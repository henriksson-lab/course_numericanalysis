// C7: shadow prices by re-solving, and a scan of the oxygen limit (minilp).
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

fn growth(glc: f64, o2: f64) -> f64 {
    let (lo, hi) = bounds(glc, o2);
    fba(lo, hi, BIO, Maximize).unwrap().0
}

fn main() {
    let h = 0.01;
    let price = |glc: f64, o2: f64| {
        let g0 = growth(glc, o2);
        ((growth(glc + h, o2) - g0) / h, (growth(glc, o2 + h) - g0) / h)
    };
    let (pg, po) = price(10.0, 15.0);
    println!("{:.3} {:.3}", pg, po + 0.0); // glucose, oxygen
    let (pg, po) = price(10.0, 4.0);
    println!("{:.3} {:.3}", pg, po);
    let (lo, hi) = bounds(10.0, 4.0);
    show(&fba(lo, hi, BIO, Maximize).unwrap().1);
    println!(" o2  growth  acetate  price_o2");
    for u in (0..15).step_by(2) {
        let (lo, hi) = bounds(10.0, u as f64);
        let (g, v) = fba(lo, hi, BIO, Maximize).unwrap();
        println!("{:3} {:7.3} {:8.3} {:9.3}", u, g, v[7] + 0.0, price(10.0, u as f64).1 + 0.0);
    }
}
