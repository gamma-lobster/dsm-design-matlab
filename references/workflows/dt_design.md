# DT DSM Design

Use this guide when the task is a standard discrete-time delta-sigma modulator design in MATLAB.

Use the native toolbox flow:

```matlab
synthesizeNTF -> realizeNTF -> stuffABCD -> simulate/analyze
```

## Core Sequence

1. Choose `order`, `OSR`, `H_inf`, `f0`, `opt`, topology, sample rate, and quantizer resolution.
2. Synthesize the NTF with `synthesizeNTF(order, OSR, opt, H_inf, f0)`.
3. Realize coefficients with `realizeNTF(ntf, form)`.
4. Build the state-space model with `stuffABCD(a, g, b, c, form)`.
5. Simulate the loop and quantizer behavior with either `simulateDSM(...)` or the explicit MATLAB loops used in the bundled examples.
6. Compute in-band SNR or SNDR from the output spectrum.
7. Check stability, state swing, and whether the achieved performance matches the spec.

Keep designs normalized unless the user gives physical circuit scaling requirements.

## Default Heuristics

| Parameter | Typical default | Notes |
| --- | --- | --- |
| `order` | 2 to 4 | Start lower if stability is uncertain |
| `OSR` | 32 or 64 | Higher OSR improves in-band noise |
| `H_inf` | 1.5 to 2.0 | Keep conservative for lowpass stability |
| `f0` | 0 | Use `0.25` for bandpass-at-`fs/4` cases |
| `opt` | 1 | Prefer optimized zeros |
| `form` | `CIFF` | Good default for lowpass multi-bit examples |
| `n_bits` | 1 to 5 | Multi-bit improves SNR but makes DAC mismatch relevant |

## Topology Guidance

- `CIFF`: default lowpass choice with straightforward signal transfer behavior.
- `CIFB`: use when feedback-form realization is preferred.
- `CRFB` or `CRFF`: use for resonator and bandpass-oriented designs.

## Primary References

- `references/designs/dt/dsm_quick_design.m`: best starting point for a new lowpass DSM design.
- `references/designs/dt/design_4th_order_ciff.m`: complete worked example with reporting and plots.
- `references/designs/dt/dsm_4th_order_simple.m`: headless or simplified variant.
