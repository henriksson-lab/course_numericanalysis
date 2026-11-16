// C1: Gillespie's direct method for mRNA birth-death, hand-written.
use rand::{rngs::StdRng, RngExt, SeedableRng};

const K: f64 = 10.0; // transcription, per hour
const DELTA: f64 = 1.0; // mRNA decay, per hour

// one cell from m = 0: event times and counts after each event
fn bd_run(t_end: f64, rng: &mut StdRng) -> (Vec<f64>, Vec<i64>) {
    let (mut t, mut m) = (0.0, 0i64);
    let (mut tt, mut mm) = (vec![0.0], vec![0i64]);
    while t < t_end {
        let (a1, a2) = (K, DELTA * m as f64);
        let a0 = a1 + a2;
        t += -(1.0 - rng.random::<f64>()).ln() / a0; // waiting time
        m += if rng.random::<f64>() * a0 < a1 { 1 } else { -1 };
        tt.push(t);
        mm.push(m);
    }
    (tt, mm)
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let (tt, _) = bd_run(10.0, &mut rng);
    println!("{}", tt.len() - 1); // number of reactions in 10 h

    let tm = [0.5, 1.0, 2.0, 10.0];
    let n = 2000;
    let mut m_at = vec![vec![0.0; n]; tm.len()];
    for c in 0..n {
        let (tt, mm) = bd_run(10.0, &mut rng);
        for (i, &t) in tm.iter().enumerate() {
            let idx = tt.partition_point(|&s| s <= t) - 1; // last event <= t
            m_at[i][c] = mm[idx] as f64;
        }
    }
    for (i, &t) in tm.iter().enumerate() {
        let v = &m_at[i];
        let mean = v.iter().sum::<f64>() / n as f64;
        let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
        println!("t {:4.1}  mean {:.4}  se {:.4}  ode {:.4}", t, mean, (var / n as f64).sqrt(),
                 10.0 * (1.0 - (-t as f64).exp()));
    }
    let v = &m_at[3];
    let mean = v.iter().sum::<f64>() / n as f64;
    let var = v.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (n - 1) as f64;
    println!("var {:.4}  fano {:.4}", var, var / mean);
    for k in [5.0, 10.0, 15.0] {
        let p_sim = v.iter().filter(|&&x| x == k).count() as f64 / n as f64;
        // Poisson(10) probability, exp(k log 10 - 10 - log k!)
        let lfact: f64 = (1..=k as u32).map(|i| (i as f64).ln()).sum();
        let p_pois = (k * 10f64.ln() - 10.0 - lfact).exp();
        println!("P(m = {:2})  sim {:.4}  poisson {:.4}", k, p_sim, p_pois);
    }
}
