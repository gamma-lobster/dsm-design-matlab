# SNR And ENOB Analysis

Use this guide when computing FFT-based SNR, SNDR, or ENOB.

## Rules

Measure performance from the quantized output spectrum, not from an internally derived noise-only waveform.

When finding the fundamental:

- Search only inside the in-band region.
- Use the output FFT after windowing.
- Exclude DC from the search.

When calculating noise:

- Exclude the identified signal bins.
- Exclude harmonic bins if the user wants SNR.
- Include harmonic bins if the user wants SNDR.
- With Hann windowing and low-frequency offsets, consider excluding bin 2 when leakage contaminates the baseband.

Use:

```matlab
fB = fs/(2*OSR);
fB_bins = ceil(N/(2*OSR));
```

Compute ENOB from:

```matlab
ENOB = (SNR - 1.76)/6.02;
```

## Primary References

- `references/analysis/debug/debug_snr.m`: use when SNR calculations disagree with expectations.
