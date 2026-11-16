// C5: the Christmas-tree gradient generator. Outlet staircase from the mixing stages,
// then the 8 streams flowing side by side in a 400 um outlet channel.
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

// n streams -> n + 1 streams: channel k takes k/n of stream k-1 and (n-k)/n of stream k
fn stage(c: &[f64]) -> Vec<f64> {
    let n = c.len();
    (0..=n).map(|k| {
        let left = if k > 0 { c[k - 1] } else { 0.0 };
        let right = if k < n { c[k] } else { 0.0 };
        (k as f64 / n as f64) * left + ((n - k) as f64 / n as f64) * right
    }).collect()
}

fn tree(inlets: &[f64], nout: usize) -> Vec<f64> {
    let mut c = inlets.to_vec();
    while c.len() < nout { c = stage(&c); }
    c
}

fn show(v: &[f64]) -> String { v.iter().map(|x| format!("{:.3}", x)).collect::<Vec<_>>().join(" ") }

fn main() {
    let lin = tree(&[0.0, 1.0], 8);
    let hump = tree(&[0.0, 1.0, 0.0], 8);
    println!("linear: {}", show(&lin));
    println!("hump:   {}", show(&hump));
    let (wo, m) = (400.0, 160);
    let dy = wo / m as f64;
    let y: Vec<f64> = (1..=m).map(|j| (j as f64 - 0.5) * dy).collect();
    let u: Vec<f64> = y.iter().map(|yj| 6.0 * 200.0 * (yj / wo) * (1.0 - yj / wo)).collect(); // mean 200 um/s
    let idx = [9, 49, 79, 109, 149]; // y = 23.75, 123.75, ..., 373.75 um
    for cs in [&lin, &hump] {
        let mut c: Vec<f64> = y.iter().map(|yj| cs[(yj / 50.0).ceil() as usize - 1]).collect(); // 8 streams of 50 um
        let mut xprev = 0.0;
        for x in [0.0, 1000.0, 5000.0, 20000.0] {
            c = march(&c, &u, 100.0, dy, 50.0, ((x - xprev) / 50.0) as usize);
            xprev = x;
            let sel: Vec<f64> = idx.iter().map(|&i| c[i]).collect();
            let (mx, mn) = (c.iter().cloned().fold(f64::MIN, f64::max), c.iter().cloned().fold(f64::MAX, f64::min));
            println!("x = {:5} um: {}   max - min {:.3}", x, show(&sel), mx - mn);
        }
    }
}
