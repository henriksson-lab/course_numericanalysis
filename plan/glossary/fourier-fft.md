**Trigonometric interpolation** — Interpolating equispaced samples of periodic data with a sum of sines and cosines with whole numbers of cycles per record.
**Discrete Fourier transform (DFT)** — X_k = Σ_j x_j e^{-2πijk/n}: the coefficients of the trigonometric interpolant, in complex form.
**DFT matrix** — The n×n matrix F with entries w^{jk}, w = e^{-2πi/n}; its columns are orthogonal and F⁻¹ = F̄ᵀ/n.
**Fast Fourier transform (FFT)** — Algorithm computing the DFT in O(n log n) by recursively splitting even and odd samples (Cooley–Tukey).
**Twiddle factor** — The factor w^k that combines the even and odd half-size DFTs in an FFT.
**Unnormalised transform** — Convention (R's `fft`, numpy's `fft`, rustfft) with no 1/n factor in the forward transform; R's inverse also omits it.
**Frequency bin** — One DFT coefficient X_k, belonging to k cycles per record (period T/k).
**Amplitude spectrum** — The amplitudes 2|X_k|/n plotted against frequency.
**Phase** — The position of a cosine in its cycle, −arg X_k; gives the peak time of a rhythm.
**Cosinor** — Least-squares fit of a mean plus one cosine and sine of known period, used to estimate amplitude and peak time of a rhythm.
**Sampling interval** — The time Δt between successive samples.
**Nyquist frequency** — f_N = 1/(2Δt), the highest frequency equispaced samples can represent.
**Aliasing** — A frequency above Nyquist producing exactly the same samples as a lower one, so it is reported at the wrong period.
**Anti-aliasing filter** — A low-pass filter applied before sampling to remove frequencies above Nyquist.
**Frequency resolution** — The spacing 1/T between DFT frequencies, set by the record length T.
**Spectral leakage** — Spreading of a peak over all frequencies when the record does not hold a whole number of cycles.
**Window (Hann window)** — A taper, w_j = ½ − ½cos(2πj/n) for Hann, multiplied into the data before the FFT to reduce leakage.
**Zero-padding** — Appending zeros before the FFT: gives a finer frequency grid or a linear instead of circular convolution, but no extra resolution.
**Low-pass filter** — Keeping only low frequencies (slow variation) and removing high ones (noise).
**Gibbs phenomenon (ringing)** — Oscillation near a jump when it is represented by a truncated sum of sines and cosines.
**Convolution** — (h * y)_j = Σ_m h_m y_{j−m}; a moving average, a blur.
**Circular convolution** — Convolution in which indices wrap around the end of the record; what a product of DFTs gives.
**Convolution theorem** — The DFT of a (circular) convolution is the entrywise product of the DFTs.
**Gain (transfer function)** — |H_k|, the factor by which a filter scales the amplitude at frequency k.
**Point spread function (PSF)** — The image of a single point source; a microscope image is the true image convolved with it.
**Deconvolution** — Undoing a convolution by dividing by the transform of the kernel; ill-conditioned when that transform is small.
**Bluestein's algorithm** — An FFT method for any length, including primes, in O(n log n).

## Install

- Nothing beyond base R (`stats::fft`, `mvfft`, `nextn`, `filter`, `convolve`), numpy/scipy (`numpy.fft`, `scipy.fft.next_fast_len`) and the Rust crates already listed (`rustfft`, `nalgebra`, `rand`, `rand_distr`).
