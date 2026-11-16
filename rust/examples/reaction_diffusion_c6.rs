// C6: dispersion relation of the Schnakenberg model; Turing threshold d_c.

fn main() {
    let (a, b) = (0.1_f64, 0.9_f64);
    let (us, vs) = (a + b, b / ((a + b) * (a + b))); // uniform steady state (1, 0.9)
    let (fu, fv, gu, gv) = (-1.0 + 2.0 * us * vs, us * us, -2.0 * us * vs, -us * us);
    // J_q = J - q^2 diag(1, d); 2 x 2, so use trace and determinant
    let tr = |q: f64, d: f64| fu + gv - q * q * (1.0 + d);
    let det = |q: f64, d: f64| (fu - q * q) * (gv - d * q * q) - fv * gu;
    let growth = |q: f64, d: f64| {
        let (t, dt) = (tr(q, d), det(q, d));
        let disc = t * t - 4.0 * dt;
        if disc < 0.0 { t / 2.0 } else { (t + disc.sqrt()) / 2.0 } // largest Re(lambda)
    };
    let q: Vec<f64> = (0..=1500).map(|i| i as f64 * 0.001).collect();
    for d in [5.0, 10.0, 20.0] {
        let g: Vec<f64> = q.iter().map(|&qq| growth(qq, d)).collect();
        let (im, gmax) = g.iter().enumerate().fold((0, f64::MIN), |m, (i, &v)| if v > m.1 { (i, v) } else { m });
        if gmax <= 0.0 {
            println!("d = {:2}: stable at every q (max growth {:.3})", d, gmax);
            continue;
        }
        let pos: Vec<f64> = q.iter().zip(&g).filter(|(_, v)| **v > 0.0).map(|(qq, _)| *qq).collect();
        println!("d = {:2}: unstable for q in [{:.3}, {:.3}], fastest q = {:.3} (wavelength {:.2}), rate {:.4}",
                 d, pos[0], pos[pos.len() - 1], q[im], 2.0 * std::f64::consts::PI / q[im], gmax);
    }
    // Turing threshold: the smallest det over q touches 0. Bisection on d.
    let min_det = |d: f64| q.iter().map(|&qq| det(qq, d)).fold(f64::MAX, f64::min);
    let (mut lo, mut hi) = (5.0, 20.0);
    while hi - lo > 1e-10 {
        let mid = 0.5 * (lo + hi);
        if min_det(mid) > 0.0 { lo = mid } else { hi = mid }
    }
    println!("{}", 0.5 * (lo + hi));
}
