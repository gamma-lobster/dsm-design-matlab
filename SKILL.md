---
name: dsm-design-matlab
description: Design and analyze delta-sigma modulators in MATLAB using Richard Schreier's Delta Sigma Toolbox. Use when Codex needs to synthesize an NTF, realize a DSM topology, build an ABCD matrix, simulate SNR or ENOB, troubleshoot stability, or model multi-bit flash ADC and thermometer DAC behavior for DSM work.
---

# DSM Design In MATLAB

Use the native MATLAB Delta Sigma Toolbox workflow:

```matlab
synthesizeNTF -> realizeNTF -> stuffABCD -> simulate/analyze
```

Prefer the bundled MATLAB references in `references/` over inventing new implementations.

## Start Here

Confirm the toolbox is on the MATLAB path before doing design work.

```matlab
% If the script is in repo root:
addpath(fullfile(fileparts(mfilename('fullpath')), 'references', 'dstoolbox'));

% If the script is inside references/:
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));

which synthesizeNTF
which realizeNTF
which stuffABCD
which simulateDSM
```

Match the `addpath(...)` pattern to the script location so the repo stays portable.

## Core Workflow

Follow this sequence unless the user asks for a narrower task.

1. Choose `order`, `OSR`, `H_inf`, `f0`, `opt`, topology, sample rate, and quantizer resolution.
2. Synthesize the NTF with `synthesizeNTF(order, OSR, opt, H_inf, f0)`.
3. Realize coefficients with `realizeNTF(ntf, form)`.
4. Build the state-space model with `stuffABCD(a, g, b, c, form)`.
5. Simulate the loop and quantizer behavior with either `simulateDSM(...)` or the explicit MATLAB loops used in the bundled examples.
6. Compute in-band SNR or SNDR from the output spectrum.
7. Check stability, state swing, and whether the achieved performance matches the spec.

Keep designs normalized unless the user gives physical circuit scaling requirements.

## Parameter Heuristics

Use these as defaults when the user does not provide a spec:

| Parameter | Typical default | Notes |
| --- | --- | --- |
| `order` | 2 to 4 | Start lower if stability is uncertain |
| `OSR` | 32 or 64 | Higher OSR improves in-band noise |
| `H_inf` | 1.5 to 2.0 | Keep conservative for lowpass stability |
| `f0` | 0 | Use `0.25` for bandpass-at-`fs/4` cases |
| `opt` | 1 | Prefer optimized zeros |
| `form` | `CIFF` | Good default for lowpass multi-bit examples |
| `n_bits` | 1 to 5 | Multi-bit improves SNR but makes DAC mismatch relevant |

Topology guidance:

- `CIFF`: default lowpass choice with straightforward signal transfer behavior.
- `CIFB`: use when feedback-form realization is preferred.
- `CRFB` or `CRFF`: use for resonator and bandpass-oriented designs.

## Simulation And Analysis Rules

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

## Stability And Debugging

If a DSM diverges or underperforms, check these first:

- Input amplitude is too high for the chosen order and `H_inf`.
- `H_inf` is too aggressive for the architecture.
- The chosen topology does not match the intended lowpass or bandpass behavior.
- The signal bin is outside the intended passband.
- The FFT and window normalization are inconsistent.
- Harmonics or DC leakage are being counted as noise incorrectly.

When debugging, report:

- coefficient vectors `a`, `g`, `b`, `c`
- `ABCD`
- peak internal state magnitudes
- output range and clipping behavior
- in-band signal bins and excluded noise bins

## Use The Bundled References

Open only the reference files needed for the current task.

- `references/dsm_quick_design.m`: best starting point for a new lowpass DSM design.
- `references/design_4th_order_ciff.m`: complete worked example with reporting and plots.
- `references/dsm_4th_order_simple.m`: headless or simplified variant.
- `references/build_dsm_simulink_model.m`: use when the user wants a Simulink model assembled from the design.
- `references/debug_snr.m`: use when SNR calculations disagree with expectations.
- `references/flash_adc_quantizer.m`: flash ADC with thermometer-code output.
- `references/thermometer_dac.m`: ideal thermometer DAC reconstruction.
- `references/thermometer_dac_mismatch.m`: static unit-cell mismatch model.
- `references/thermometer_dac_dwa.m`: dynamic weighted averaging for mismatch shaping.
- `references/design_4th_order_ciff_with_dac_mismatch.m`: end-to-end mismatch and DWA example.

If the user asks about toolbox behavior, inspect the relevant source in `references/dstoolbox/` rather than guessing.

## Component-Level Modeling

Support these adjacent tasks when they are part of the DSM design problem:

- model a flash ADC that outputs thermometer code
- reconstruct the analog level with a thermometer DAC
- inject unit-cell mismatch into a multi-bit DAC
- compare ideal, mismatched, and DWA-shaped feedback DAC behavior

Use the bundled component functions before writing new helpers. Keep thermometer outputs contiguous unless the task explicitly requires rotation or scrambling. For DWA, rotate cell usage for DAC feedback while preserving the ADC's original code meaning.

## Response Pattern

When helping with a design, prefer this output structure:

1. Restate the target spec in MATLAB terms.
2. Show the chosen design parameters and why they are reasonable.
3. Provide or edit MATLAB code that follows the toolbox workflow.
4. Summarize the expected stability and performance limits.
5. Call out any assumptions, especially around OSR, topology, quantizer bits, and mismatch.

Keep explanations practical and tied to the provided MATLAB files.
