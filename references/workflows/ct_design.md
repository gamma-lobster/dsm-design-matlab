# CT DSM Design

Use this guide for continuous-time delta-sigma modulator design, mapping, and simulation.

Use the CT workflow:

```matlab
synthesizeNTF -> realizeNTF_ct -> mapCtoD -> simulate/analyze
```

## Core Sequence

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

## Default Heuristics

| Parameter | Typical default | Notes |
| --- | --- | --- |
| `ct_form` | `FF` | Good default when mapping a lowpass CIFF-like design to CT |
| `tdac` | `[0 1]` | Use `[0.5 1.5]` to study half-cycle ELD with NRZ feedback |
| `jitter_mode` | `edge` | Use `equivalent` for the sampled textbook-style DAC jitter model |

## CT Guidance

- CT `FF`: good match when you want a CIFF-like CT realization with feedforward output summation.
- CT `FB`: use when a feedback-form CT loop filter is specifically desired.
- Delayed `tdac`: use to inject excess loop delay directly into the CT realization and Simulink DAC waveform.
- `jitter_mode='edge'`: continuous-time DAC edge perturbation, useful for waveform-level jitter studies.
- `jitter_mode='equivalent'`: sampled DAC jitter error injection, useful when comparing against equation-based SJNR estimates.

## Primary References

- `references/designs/ct/design_3rd_order_ct_dsm_10mhz.m`: best starting point for the current continuous-time 3rd-order reference flow.
- `references/simulink/builders/build_ct_dsm_simulink_model.m`: generic CTDSM Simulink model with CT integrators and a sampled quantizer input.
- `references/simulink/runs/run_3rd_order_ct_simulink_model.m`: execute and plot the CT Simulink model.
- `references/simulink/runs/run_3rd_order_ct_eld_simulink_model.m`: verified uncompensated-versus-compensated ELD comparison.
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_simulink_model.m`: DAC clock-jitter sweep in the CT Simulink model.
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_theory_compare.m`: compare edge-jitter, equivalent-jitter, and equation-style SJNR side by side.
