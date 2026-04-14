# Stability And Debugging

Use this guide when a DSM diverges, clips, or underperforms.

## Check These First

- Input amplitude is too high for the chosen order and `H_inf`.
- `H_inf` is too aggressive for the architecture.
- The chosen topology does not match the intended lowpass or bandpass behavior.
- For CTDSM, the `ct_form` does not match the intended CIFF-like or CIFB-like structure.
- For CTDSM, the sampler is misplaced or missing ahead of the quantizer.
- For CTDSM, the integrator drives are missing the required `fs` scaling.
- For CTDSM, the DAC timing used in Simulink does not match the `tdac` used in `realizeNTF_ct`.
- For CTDSM with ELD, compensation may be missing if delayed `tdac` is applied to an uncompensated loop filter.
- For CTDSM jitter studies, confirm whether the run is using edge-time jitter or equivalent sampled jitter before comparing to theory.
- The signal bin is outside the intended passband.
- The FFT and window normalization are inconsistent.
- In coherent FFT tests, remember a tone at DFT bin `k` appears at MATLAB FFT array index `k+1`.
- Harmonics or DC leakage are being counted as noise incorrectly.

## Report These When Debugging

- coefficient vectors `a`, `g`, `b`, `c`
- `ABCD`
- `ABCDc`, `tdac`, and `tdac2` for CT realizations
- whether the CT case is uncompensated ELD or compensated ELD
- peak internal state magnitudes
- output range and clipping behavior
- in-band signal bins and excluded noise bins
