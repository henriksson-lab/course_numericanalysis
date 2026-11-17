// C1: diffusion coefficient of a bead from the mean squared displacement (MSD) of its track.
use rand::{rngs::StdRng, SeedableRng};
use rand_distr::{Distribution, Normal};

const D: f64 = 0.43; // µm^2/s
const DT: f64 = 0.1; // s per frame

// n frames of a 2D track (columns x, y), flow u along x in µm/s
fn track(n: usize, u: f64, rng: &mut StdRng) -> Vec<[f64; 2]> {
    let step = Normal::new(0.0, (2.0 * D * DT).sqrt()).unwrap();
    let (mut x, mut y) = (0.0, 0.0);
    (1..=n)
        .map(|i| {
            x += step.sample(rng);
            y += step.sample(rng);
            [x + u * DT * i as f64, y]
        })
        .collect()
}

// D from a line through 0 fitted to MSD at lags 1..K; dims = 2 (x, y) or 1 (y only)
fn dhat(xy: &[[f64; 2]], k_max: usize, dims: usize) -> f64 {
    let (mut stm, mut stt) = (0.0, 0.0);
    for k in 1..=k_max {
        let n = xy.len() - k;
        let mut msd = 0.0;
        for i in 0..n {
            let dy = xy[i + k][1] - xy[i][1];
            let dx = if dims == 2 { xy[i + k][0] - xy[i][0] } else { 0.0 };
            msd += dx * dx + dy * dy;
        }
        msd /= n as f64;
        let lag = k as f64 * DT;
        stm += lag * msd;
        stt += lag * lag;
    }
    stm / stt / (2.0 * dims as f64) // slope / (2 d)
}

fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len() as f64;
    let m = v.iter().sum::<f64>() / n;
    (m, (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (n - 1.0)).sqrt())
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    let xy = track(100, 0.0, &mut rng);
    println!("K4 {:.4}  K50 {:.4}", dhat(&xy, 4, 2), dhat(&xy, 50, 2));
    for (n, k) in [(100, 4), (100, 50), (1000, 4)] {
        let d: Vec<f64> = (0..1000).map(|_| dhat(&track(n, 0.0, &mut rng), k, 2)).collect();
        let (m, s) = mean_sd(&d);
        println!("n = {:4}  K = {:2}  mean {:.3}  sd {:.3}", n, k, m, s);
    }
    let (mut d2, mut dy) = (Vec::new(), Vec::new());
    for _ in 0..1000 {
        let xy = track(100, 5.0, &mut rng);
        d2.push(dhat(&xy, 4, 2));
        dy.push(dhat(&xy, 4, 1));
    }
    println!("2D with flow {:.4}  y only {:.4}", mean_sd(&d2).0, mean_sd(&dy).0);
}
