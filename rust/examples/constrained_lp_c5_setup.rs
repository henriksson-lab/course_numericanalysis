// C5 setup: the stoichiometric matrix of the toy network (rows = metabolites).
const RXN: [&str; 9] = ["EX_glc", "EX_o2", "GLY", "BYP1", "BYP2", "RES", "FER", "EX_ac", "BIO"];
const S: [[f64; 9]; 6] = [
    [1., 0., -1., -1., 0., 0., 0., 0., 0.],   // G
    [0., 1., 0., 0., 0., -2., 0., 0., 0.],    // O2
    [0., 0., 0., 1., -1., 0., 0., 0., 0.],    // X
    [0., 0., 2., 0., 2., -1., -1., 0., -1.],  // P
    [0., 0., 2., 0., 2., 11., 1., 0., -5.],   // ATP
    [0., 0., 0., 0., 0., 0., 1., -1., 0.],    // Ac
];

fn main() {
    println!("{} x {}  ({} reactions)", S.len(), S[0].len(), RXN.len());
}
