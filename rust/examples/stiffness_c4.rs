// C4: full mass-action model vs the Michaelis-Menten (QSSA) model for three enzyme levels.
// Reference solutions with Dopri5 at tight tolerance (stiffness test off: slow but correct).
use ode_solvers::continuous_output_model::ContinuousOutputModel;
use ode_solvers::{Dopri5, OutputType, System, Vector1, Vector4};

const KF: f64 = 100.0; const KR: f64 = 100.0; const KCAT: f64 = 1.0;
const KM: f64 = (KR + KCAT) / KF;

struct Enzyme;
impl System<f64, Vector4<f64>> for Enzyme {
    fn system(&self, _t: f64, y: &Vector4<f64>, dy: &mut Vector4<f64>) {
        let b = KF * y[1] * y[0] - KR * y[2];
        dy[0] = -b; dy[1] = -b + KCAT * y[2]; dy[2] = b - KCAT * y[2]; dy[3] = KCAT * y[2];
    }
}
struct Qssa { e0: f64 }
impl System<f64, Vector1<f64>> for Qssa {
    fn system(&self, _t: f64, y: &Vector1<f64>, dy: &mut Vector1<f64>) {
        dy[0] = -KCAT * self.e0 * y[0] / (KM + y[0]);
    }
}

fn main() {
    for e0 in [0.1, 1.0, 5.0] {
        let tend = 20.0 / e0;                             // long enough to use up S0 = 10
        let mut full = ContinuousOutputModel::default();
        Dopri5::from_param(Enzyme, 0.0, tend, tend, Vector4::new(10.0, e0, 0.0, 0.0), 1e-10, 1e-12,
                           0.9, 0.04, 0.2, 10.0, tend, 0.0, 10_000_000, u32::MAX, OutputType::Dense)
            .integrate_with_continuous_output_model(&mut full).unwrap();
        let mut qssa = ContinuousOutputModel::default();
        Dopri5::new(Qssa { e0 }, 0.0, tend, tend, Vector1::new(10.0), 1e-10, 1e-12)
            .integrate_with_continuous_output_model(&mut qssa).unwrap();
        let err = (0..=400).map(|k| {
            let t = tend * k as f64 / 400.0;
            (full.evaluate(t).unwrap()[3] - (10.0 - qssa.evaluate(t).unwrap()[0])).abs()
        }).fold(0.0, f64::max);
        println!("E0 = {:3.1}  max |P_full - P_qssa| = {:.4} uM  (E0 / (S0 + Km) = {:.3})",
                 e0, err, e0 / (10.0 + KM));
    }
}
