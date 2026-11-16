fn main() {
    let len: [i32; 13] = [248956422, 242193529, 198295559, 190214555, 181538259,
                          170805979, 159345973, 145138636, 138394717, 133797422,
                          135086622, 133275309, 114364328];   // chr1-13, GRCh38

    // Plain `total + l` would stop the program in a debug build
    // ("attempt to add with overflow") and wrap silently in a release build.
    // Rust makes you choose explicitly:
    let mut checked: Option<i32> = Some(0); // checked_add: None on overflow
    let mut wrapped: i32 = 0;               // wrapping_add: wraps around like numpy int32
    let mut as_f64: f64 = 0.0;              // convert to double: the fix
    for (chr, &l) in len.iter().enumerate() {
        checked = checked.and_then(|t| t.checked_add(l));
        wrapped = wrapped.wrapping_add(l);
        as_f64 += l as f64;
        if chr >= 10 {                      // chromosomes 11-13 (index starts at 0)
            println!("chr{}: {:?} {} {}", chr + 1, checked, wrapped, as_f64);
        }
    }
}
