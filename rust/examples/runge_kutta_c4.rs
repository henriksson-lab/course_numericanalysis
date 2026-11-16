// C4: an adaptive Heun-Euler 2(1) solver, written by hand
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
    // (a) logistic to t = 10
    let f = |_t: f64, n: &[f64]| vec![0.8 * n[0] * (1.0 - n[0] / 1.2)];
    let exact = |t: f64| 1.2 / (1.0 + (1.2 / 0.05 - 1.0) * (-0.8 * t).exp());
    for rtol in [1e-3, 1e-6] {
        let s = adapt_heun(&f, 0.0, &[0.05], 10.0, 1.0, rtol, 1e-9);
        let n = s.h.len();
        println!("rtol {:e}: error {:.3e}, steps {}, rejected {}, f evaluations {}",
                 rtol, s.y[n][0] - exact(10.0), n, s.nrej, 2 * (n + s.nrej));
    }
    // (b) SIR epidemic to day 200
    let (beta, gamma) = (0.3, 0.1);
    let sir = |_t: f64, y: &[f64]| vec![-beta * y[0] * y[1], beta * y[0] * y[1] - gamma * y[1], gamma * y[1]];
    let s = adapt_heun(&sir, 0.0, &[0.999, 0.001, 0.0], 200.0, 0.1, 1e-6, 1e-9);
    println!("steps {}, rejected {}", s.h.len(), s.nrej);
    for tq in [0.0, 13.0, 38.0, 100.0, 199.0] {
        let i = s.t.iter().filter(|&&ti| ti <= tq).count() - 1; // step that starts before tq
        println!("h near t = {}: {:.3e}", tq, s.h[i]);
    }
}
