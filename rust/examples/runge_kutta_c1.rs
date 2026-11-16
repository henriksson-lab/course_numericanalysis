// C1: midpoint, Heun and RK4 by hand; order on the logistic curve
type Rhs = dyn Fn(f64, &[f64]) -> Vec<f64>;
type Step = fn(&Rhs, f64, &[f64], f64) -> Vec<f64>;

// y + a * k, element by element (R and numpy do this with `y + a * k`)
fn axpy(y: &[f64], a: f64, k: &[f64]) -> Vec<f64> {
    y.iter().zip(k).map(|(yi, ki)| yi + a * ki).collect()
}

fn euler_step(f: &Rhs, t: f64, y: &[f64], h: f64) -> Vec<f64> {
    axpy(y, h, &f(t, y))
}
fn midpoint_step(f: &Rhs, t: f64, y: &[f64], h: f64) -> Vec<f64> {
    let k1 = f(t, y);
    let k2 = f(t + h / 2.0, &axpy(y, h / 2.0, &k1));
    axpy(y, h, &k2)
}
fn heun_step(f: &Rhs, t: f64, y: &[f64], h: f64) -> Vec<f64> {
    let k1 = f(t, y);
    let k2 = f(t + h, &axpy(y, h, &k1));
    let k: Vec<f64> = k1.iter().zip(&k2).map(|(a, b)| a + b).collect();
    axpy(y, h / 2.0, &k)
}
fn rk4_step(f: &Rhs, t: f64, y: &[f64], h: f64) -> Vec<f64> {
    let k1 = f(t, y);
    let k2 = f(t + h / 2.0, &axpy(y, h / 2.0, &k1));
    let k3 = f(t + h / 2.0, &axpy(y, h / 2.0, &k2));
    let k4 = f(t + h, &axpy(y, h, &k3));
    let k: Vec<f64> = (0..y.len()).map(|i| k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]).collect();
    axpy(y, h / 6.0, &k)
}
// n fixed steps of any method; returns the final state
fn rk_solve(step: Step, f: &Rhs, t0: f64, y0: &[f64], h: f64, n: usize) -> Vec<f64> {
    let mut y = y0.to_vec();
    for i in 0..n {
        y = step(f, t0 + i as f64 * h, &y, h);
    }
    y
}

fn main() {
    let methods: [(&str, Step); 4] = [("Euler", euler_step), ("Midpoint", midpoint_step),
                                      ("Heun", heun_step), ("RK4", rk4_step)];
    // (a) check B1 and B2: drug clearance C' = -0.8 C, C0 = 10 mg/L, one step h = 0.5
    let clear = |_t: f64, c: &[f64]| vec![-0.8 * c[0]];
    for (name, s) in methods {
        println!("{}: {:.7}", name, s(&clear, 0.0, &[10.0], 0.5)[0]);
    }
    println!("exact: {:.7}", 10.0 * (-0.4f64).exp());

    // (b) logistic: error in N(10) for h = 1, 1/2, ..., 1/64
    let f = |_t: f64, n: &[f64]| vec![0.8 * n[0] * (1.0 - n[0] / 1.2)];
    let exact = |t: f64| 1.2 / (1.0 + (1.2 / 0.05 - 1.0) * (-0.8 * t).exp());
    let mut err = vec![[0.0; 3]; 7];
    for (i, row) in err.iter_mut().enumerate() {
        let h = 0.5f64.powi(i as i32);
        let n = (10.0 / h) as usize;
        for (j, s) in [euler_step as Step, heun_step, rk4_step].into_iter().enumerate() {
            row[j] = (rk_solve(s, &f, 0.0, &[0.05], h, n)[0] - exact(10.0)).abs();
        }
        println!("h = {:<9} {:.2e} {:.2e} {:.2e}", h, row[0], row[1], row[2]);
    }
    println!("observed order (Euler, Heun, RK4):");
    for i in 0..6 {
        let p: Vec<String> = (0..3).map(|j| format!("{:.2}", (err[i][j] / err[i + 1][j]).log2())).collect();
        println!("{}", p.join(" "));
    }
}
