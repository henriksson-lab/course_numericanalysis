fn main() {
    let k = 10.0;
    let dose = |r: f64| k * r / (1.0 - r); // closure, captures k
    for r0 in [0.5, 0.9, 0.99] {
        let lo = dose(r0 * 0.995); // response off by -0.5 %
        let hi = dose(r0 * 1.005); // and by +0.5 %
        println!("R = {:.2}   dose from {:.1} to {:.1}   ratio {:.2}", r0, lo, hi, hi / lo);
    }
}
