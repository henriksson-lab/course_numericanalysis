// C4: the SIR epidemic with a hand-written vector Euler method.
const BETA: f64 = 0.3; // per day
const GAMMA: f64 = 0.1; // per day
const N: f64 = 10000.0; // population size

fn sir(_t: f64, y: &[f64]) -> Vec<f64> {
    let (s, i) = (y[0], y[1]);
    vec![-BETA * s * i / N, BETA * s * i / N - GAMMA * i, GAMMA * i]
}

// one Vec per time point; f takes the state as a slice &[f64]
fn euler_sys(f: impl Fn(f64, &[f64]) -> Vec<f64>, t0: f64, y0: &[f64], h: f64, n: usize) -> Vec<Vec<f64>> {
    let mut ys = vec![y0.to_vec()];
    for i in 0..n {
        let y = &ys[i];
        let slope = f(t0 + i as f64 * h, y);
        let next: Vec<f64> = y.iter().zip(&slope).map(|(yj, fj)| yj + h * fj).collect();
        ys.push(next);
    }
    ys
}

// index and value of the largest I
fn peak(ys: &[Vec<f64>]) -> (usize, f64) {
    ys.iter().enumerate().map(|(i, y)| (i, y[1])).fold((0, f64::MIN), |a, b| if b.1 > a.1 { b } else { a })
}

fn main() {
    let y0 = [9990.0, 10.0, 0.0];

    // (a) one run with h = 1 day
    let ys = euler_sys(sir, 0.0, &y0, 1.0, 160);
    let (day, imax) = peak(&ys);
    println!("peak {:.3} on day {}", imax, day);
    let drift = ys.iter().map(|y| (y[0] + y[1] + y[2] - N).abs()).fold(0.0, f64::max);
    println!("largest |S + I + R - N|: {:.1e}", drift);
    // No plotting here: write the rows to a file and plot in R or Python.

    // (b) error in the peak for several step sizes
    let r0 = BETA / GAMMA;
    let exact = y0[0] + y0[1] - N / r0 * (1.0 + (r0 * y0[0] / N).ln());
    println!("exact peak {:.3}", exact);
    for h in [4.0, 2.0, 1.0, 0.5, 0.25, 0.1, 0.01] {
        let ys = euler_sys(sir, 0.0, &y0, h, (160.0 / h).round() as usize);
        let (i, imax) = peak(&ys);
        println!("{:5.2}  {:9.3}  {:8.3}  {:6.2}", h, imax, imax - exact, i as f64 * h);
    }
}
