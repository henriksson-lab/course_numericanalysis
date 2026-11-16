// C6: early extinction of an epidemic (SIR, N = 10 000, beta = 0.3, gamma = 0.1 per day).
use rand::{rngs::StdRng, RngExt, SeedableRng};

fn sir_fate(i0: i64, rng: &mut StdRng) -> bool {
    let (n, beta, gamma) = (10000.0, 0.3, 0.1);
    let (mut s, mut i) = (10000 - i0, i0);
    while i > 0 && i < 100 { // stop: extinct, or 100 infected
        let a1 = beta * s as f64 * i as f64 / n; // infection
        let a2 = gamma * i as f64; // recovery
        if rng.random::<f64>() * (a1 + a2) < a1 { s -= 1; i += 1; } else { i -= 1; }
    }
    i == 0 // true if the outbreak died out
}

fn main() {
    let mut rng = StdRng::seed_from_u64(1);
    for i0 in 1..=3 {
        let ext = (0..10000).filter(|_| sir_fate(i0, &mut rng)).count() as f64 / 10000.0;
        println!("I0 = {}  extinct {:.4}  se {:.4}  branching (1/3)^I0 = {:.4}", i0, ext,
                 (ext * (1.0 - ext) / 9999.0).sqrt(), (1.0f64 / 3.0).powi(i0 as i32));
    }
}
