// C6: atol too large: the hand-written adaptive solver of C4 on the SIR epidemic
struct Run { t: Vec<f64>, h: Vec<f64>, y: Vec<Vec<f64>>, nrej: usize }

fn adapt_heun(f: &dyn Fn(f64, &[f64]) -> Vec<f64>, mut t: f64, y0: &[f64], tend: f64,
              mut h: f64, rtol: f64, atol: f64) -> Run {
    let mut y = y0.to_vec();
    let mut run = Run { t: vec![t], h: vec![], y: vec![y.clone()], nrej: 0 };
    let m = y.len();
    while t < tend {
        h = h.min(tend - t);
        let k1 = f(t, &y);
        let ye: Vec<f64> = (0..m).map(|i| y[i] + h * k1[i]).collect(); // Euler, order 1
        let k2 = f(t + h, &ye);
        let yhi: Vec<f64> = (0..m).map(|i| y[i] + h / 2.0 * (k1[i] + k2[i])).collect(); // Heun
        let mut ss = 0.0;
        for i in 0..m {
            let sc = atol + rtol * y[i].abs().max(yhi[i].abs());
            ss += ((yhi[i] - ye[i]) / sc).powi(2);
        }
        let err = (ss / m as f64).sqrt();
        if err <= 1.0 { // accept and record
            t += h;
            y = yhi;
            run.t.push(t); run.h.push(h); run.y.push(y.clone());
        } else { // reject: retry with smaller h
            run.nrej += 1;
        }
        h *= (0.9 * err.powf(-0.5)).clamp(0.2, 5.0);
    }
    run
}

fn main() {
    let (beta, gamma, n) = (0.3, 0.1, 10000.0);
    let sir = |_t: f64, y: &[f64]| {
        let (s, i) = (y[0], y[1]);
        vec![-beta * s * i / n, beta * s * i / n - gamma * i, gamma * i]
    };
    let i_ref = 998.2012187; // I(60), from a very tight solve
    for atol in [10.0, 1.0, 1e-6] {
        print!("atol {:e}:", atol);
        for k in 3..7 {
            let rtol = 10f64.powi(-k);
            let s = adapt_heun(&sir, 0.0, &[9990.0, 10.0, 0.0], 60.0, 0.1, rtol, atol);
            print!("  {:.1e}", (s.y.last().unwrap()[1] / i_ref - 1.0).abs());
        }
        println!();
    }
}
