// C7: AUC from the blood samples: linear trapezoid vs linear-up/log-down,
// per interval and in total, and extrapolation to infinity.
fn model(t: f64) -> f64 {
    500.0 / 50.0 * 1.5 / 1.3 * ((-0.2 * t).exp() - (-1.5 * t).exp())
}

fn main() {
    let t: [f64; 10] = [0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0];
    let c: [f64; 10] = [0.0, 4.99, 6.87, 7.16, 6.20, 5.16, 3.47, 2.33, 1.05, 0.09];
    let (mut lin_sum, mut lul_sum, mut true_sum) = (0.0, 0.0, 0.0);
    let (mut e_lin, mut e_lul) = (Vec::new(), Vec::new());
    for i in 0..9 {
        let (dt, y1, y2) = (t[i + 1] - t[i], c[i], c[i + 1]);
        let lin = dt * (y1 + y2) / 2.0;
        let lul = if y2 < y1 && y2 > 0.0 { dt * (y1 - y2) / (y1 / y2).ln() } else { lin };
        let tr = quadrature::integrate(model, t[i], t[i + 1], 1e-12).integral;
        e_lin.push(format!("{:.3}", lin - tr));
        e_lul.push(format!("{:.3}", lul - tr));
        lin_sum += lin; lul_sum += lul; true_sum += tr;
    }
    println!("lin {}", e_lin.join(" "));
    println!("lul {}", e_lul.join(" "));
    println!("{} {} {}", lin_sum, lul_sum, true_sum);
    // terminal slope: least-squares line through log(c) at 8, 12, 24 h
    let (x, y): (Vec<f64>, Vec<f64>) = (7..10).map(|i| (t[i], c[i].ln())).unzip();
    let (mx, my) = (x.iter().sum::<f64>() / 3.0, y.iter().sum::<f64>() / 3.0);
    let sxy: f64 = x.iter().zip(&y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    let lz = -sxy / sxx;
    println!("{} {}", lin_sum + c[9] / lz, lul_sum + c[9] / lz);
    // model route: no infinite limits in the crate; substitute t = s / (1 - s)
    let g = |s: f64| model(s / (1.0 - s)) / (1.0 - s).powi(2);
    println!("{}", quadrature::integrate(g, 0.0, 1.0, 1e-10).integral);
}
