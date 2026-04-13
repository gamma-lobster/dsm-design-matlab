---
name: dsm-design-matlab
description: Design and analyze delta-sigma modulators in MATLAB using Richard Schreier's Delta Sigma Toolbox. Use when Codex needs to synthesize an NTF, realize DT or CT DSM topologies, build ABCD or ABCDc matrices, assemble MASH structures, simulate SNR or ENOB, troubleshoot stability, or model multi-bit flash ADC and thermometer DAC behavior for DSM work.
---

# DSM Design In MATLAB

Use the native MATLAB Delta Sigma Toolbox workflow:

```matlab
synthesizeNTF -> realizeNTF -> stuffABCD -> simulate/analyze
```

For continuous-time DSM work, use:

```matlab
synthesizeNTF -> realizeNTF_ct -> mapCtoD -> simulate/analyze
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
which realizeNTF_ct
which stuffABCD
which mapCtoD
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

For DT MASH work:

1. Design or provide one NTF per stage.
2. Realize each stage with `realizeNTF(..., 'CIFF')` or the requested form.
3. Drive stage 1 from the external input and each later stage from the previous stage residue.
4. Form the residue as `y_k - v_k`.
5. For textbook cascaded residue-shaping MASH cancellation, combine the stage outputs digitally with cumulative previous-stage NTF products:
   `v_mash = v_1 + NTF_1(z)v_2 + NTF_1(z)NTF_2(z)v_3 + ...`
6. Prefer explicit stage-by-stage simulation when debugging residue transfer or cancellation.

For CTDSM work:

1. Synthesize the target DT NTF with `synthesizeNTF(...)`.
2. Realize the CT loop filter with `realizeNTF_ct(ntf, ct_form, tdac)`.
3. Use `ct_form='FF'` when a CIFF-like feedforward CT architecture is desired.
4. Convert the CT loop to a sampled equivalent with `mapCtoD(...)` for sample-time validation.
5. When building a CT Simulink model with continuous integrators, scale each integrator drive by `fs` (equivalently `1/Ts`) so the normalized toolbox coefficients map correctly into seconds.
6. Place an explicit sampler before the ADC or quantizer.
7. Model excess loop delay by adjusting `tdac`, for example `tdac=[0.5 1.5]` for half-cycle ELD with an NRZ DAC.
8. For compensated ELD cases, let `realizeNTF_ct(...)` synthesize the needed direct or compensation path from the delayed `tdac`.
9. For DAC clock jitter on a delayed NRZ-like CT feedback pulse, start with `jitter_mode='edge'` to perturb the DAC update edge directly in Simulink.
10. If you want a sampled textbook-style DAC jitter model, use `jitter_mode='equivalent'`, which injects an error term proportional to `v[n] - v[n-1]`.
11. Use `references/analysis/jitter/calculate_ct_dac_jitter_sjnr.m` for the current standalone equation-style `J` and `SJNR` calculation, including the `2*pi` conversion from the numerically evaluated `f`-domain integral to the `w`-domain form of Eq. 9.17.
12. Measure SNR from the sampled quantizer output, and inspect inter-sample waveforms separately.

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
| `MASH structure` | `2-1` | Good first cascaded DT reference example |
| `ct_form` | `FF` | Good default when mapping a lowpass CIFF-like design to CT |
| `tdac` | `[0 1]` | Use `[0.5 1.5]` to study half-cycle ELD with NRZ feedback |
| `jitter_mode` | `edge` | Use `equivalent` for the sampled textbook-style DAC jitter model |
| `n_bits` | 1 to 5 | Multi-bit improves SNR but makes DAC mismatch relevant |

Topology guidance:

- `CIFF`: default lowpass choice with straightforward signal transfer behavior.
- `MASH 2-1`: good introductory DT cascaded-noise-shaping case when you want a 3rd-order overall response from stable lower-order stages.
- `CIFB`: use when feedback-form realization is preferred.
- `CRFB` or `CRFF`: use for resonator and bandpass-oriented designs.
- CT `FF`: good match when you want a CIFF-like CT realization with feedforward output summation.
- CT `FB`: use when a feedback-form CT loop filter is specifically desired.
- Delayed `tdac`: use to inject excess loop delay directly into the CT realization and Simulink DAC waveform.
- `jitter_mode='edge'`: continuous-time DAC edge perturbation, useful for waveform-level jitter studies.
- `jitter_mode='equivalent'`: sampled DAC jitter error injection, useful when comparing against equation-based SJNR estimates.

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

When debugging, report:

- coefficient vectors `a`, `g`, `b`, `c`
- `ABCD`
- `ABCDc`, `tdac`, and `tdac2` for CT realizations
- whether the CT case is uncompensated ELD or compensated ELD
- peak internal state magnitudes
- output range and clipping behavior
- in-band signal bins and excluded noise bins

## Use The Bundled References

Open only the reference files needed for the current task.

- `references/designs/dt/dsm_quick_design.m`: best starting point for a new lowpass DSM design.
- `references/designs/ct/design_3rd_order_ct_dsm_10mhz.m`: best starting point for the current continuous-time 3rd-order reference flow.
- `references/designs/dt/design_4th_order_ciff.m`: complete worked example with reporting and plots.
- `references/designs/dt/design_mash_2_1_ciff_10mhz.m`: best starting point for the current DT 2-1 MASH reference flow.
- `references/designs/dt/dsm_4th_order_simple.m`: headless or simplified variant.
- `references/simulink/builders/build_dsm_simulink_model.m`: use when the user wants a Simulink model assembled from the design.
- `references/simulink/builders/build_mash_dsm_simulink_model.m`: use when the user wants a generic DT MASH Simulink model assembled from stage descriptions and digital cancellation filters.
- `references/simulink/builders/build_mash_21_dsm_simulink_model.m`: thin compatibility wrapper for the current 2-1 MASH example.
- `references/simulink/builders/build_ct_dsm_simulink_model.m`: use when the user wants a generic CTDSM Simulink model with CT integrators and a sampled quantizer input.
- `references/simulink/runs/run_mash_2_1_simulink_model.m`: use when the user wants the current MASH Simulink model executed and plotted.
- `references/simulink/runs/run_3rd_order_ct_simulink_model.m`: use when the user wants the CT Simulink model executed and plotted.
- `references/simulink/runs/run_3rd_order_ct_eld_simulink_model.m`: use when the user wants a verified uncompensated-versus-compensated ELD comparison in the CT Simulink model.
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_simulink_model.m`: use when the user wants a DAC clock-jitter sweep in the CT Simulink model.
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_theory_compare.m`: use when the user wants edge-jitter, equivalent-jitter, and equation-style SJNR compared side by side.
- `references/analysis/jitter/calculate_ct_dac_jitter_sjnr.m`: use when the user wants the standalone `J` and `SJNR` calculation.
- `references/analysis/debug/debug_snr.m`: use when SNR calculations disagree with expectations.
- `references/components/flash_adc_quantizer.m`: flash ADC with thermometer-code output.
- `references/components/thermometer_dac.m`: ideal thermometer DAC reconstruction.
- `references/components/thermometer_dac_mismatch.m`: static unit-cell mismatch model.
- `references/components/thermometer_dac_dwa.m`: dynamic weighted averaging for mismatch shaping.
- `references/designs/dt/design_4th_order_ciff_with_dac_mismatch.m`: end-to-end mismatch and DWA example.

If the user asks about toolbox behavior, inspect the relevant source in `references/third_party/dstoolbox/` rather than guessing.

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
3. Provide or edit MATLAB code that follows the DT or CT toolbox workflow.
4. Summarize the expected stability and performance limits.
5. Call out any assumptions, especially around OSR, topology, stage partitioning for MASH, quantizer bits, DAC timing, ELD compensation, and mismatch.

Keep explanations practical and tied to the provided MATLAB files.
