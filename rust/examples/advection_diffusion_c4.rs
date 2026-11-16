// C4: the H-filter. Fraction of the solute flux in the buffer half vs distance, for a
// small molecule, a protein and a cell; the shortest channel that extracts 45 % of the
// small molecule.
fn thomas(a: &[f64], b: &[f64], c: &[f64], d: &[f64]) -> Vec<f64> {
    let n = d.len();
    let (mut b, mut d) = (b.to_vec(), d.to_vec());
    for i in 1..n {
        let w = a[i] / b[i - 1];
        b[i] -= w * c[i - 1];
        d[i] -= w * d[i - 1];
    }
    let mut x = vec![0.0; n];
    x[n - 1] = d[n - 1] / b[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
    }
    x
}

// backward Euler in x (as in C3)
fn march(c: &[f64], u: &[f64], d: f64, dy: f64, dx: f64, nsteps: usize) -> Vec<f64> {
    let m = c.len();
    let a: Vec<f64> = u.iter().map(|uj| d * dx / (uj * dy * dy)).collect();
    let mut di: Vec<f64> = a.iter().map(|aj| 1.0 + 2.0 * aj).collect();
    di[0] = 1.0 + a[0];
    di[m - 1] = 1.0 + a[m - 1];
    let off: Vec<f64> = a.iter().map(|aj| -aj).collect();
    let mut c = c.to_vec();
    for _ in 0..nsteps {
        c = thomas(&off, &di, &off, &c);
    }
    c
}

fn main() {
    let (w, ubar, m) = (100.0, 1000.0, 50);
    let dy = w / m as f64;
    let y: Vec<f64> = (1..=m).map(|j| (j as f64 - 0.5) * dy).collect();
    let u: Vec<f64> = y.iter().map(|yj| 6.0 * ubar * (yj / w) * (1.0 - yj / w)).collect();
    let c0: Vec<f64> = y.iter().map(|&yj| if yj < w / 2.0 { 1.0 } else { 0.0 }).collect();
    let total: f64 = u.iter().zip(&c0).map(|(a, b)| a * b).sum();
    // fraction in the buffer half (y > w/2) at each x
    let crossed = |d: f64, xs: &[f64]| -> Vec<f64> {
        let (mut c, mut xprev, dx) = (c0.clone(), 0.0, 50.0);
        xs.iter().map(|&x| {
            c = march(&c, &u, d, dy, dx, ((x - xprev) / dx).round() as usize);
            xprev = x;
            (0..m).filter(|&j| y[j] > w / 2.0).map(|j| u[j] * c[j]).sum::<f64>() / total
        }).collect()
    };
    let xs = [1000.0, 2000.0, 5000.0, 10000.0, 20000.0];
    for d in [1000.0, 100.0, 0.05] {
        let f: Vec<String> = crossed(d, &xs).iter().map(|v| format!("{:.4}", v)).collect();
        println!("D = {:7.2}: {}", d, f.join(" "));
    }
    // shortest length (50 um steps) at which the small molecule reaches 45 %
    let xs2: Vec<f64> = (1..=100).map(|k| 50.0 * k as f64).collect();
    let (fs, fp) = (crossed(1000.0, &xs2), crossed(100.0, &xs2));
    let j = fs.iter().position(|&f| f >= 0.45).unwrap();
    println!("length {} um: small molecule {:.3}, protein {:.3}", xs2[j], fs[j], fp[j]);
}
