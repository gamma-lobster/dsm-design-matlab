# Codex_DSM

MATLAB reference workspace for delta-sigma modulator design using Richard Schreier's Delta Sigma Toolbox, plus supporting examples for flash ADC quantization, thermometer DAC feedback, DAC mismatch, DWA, and generated Simulink models.

## Current State

This repo appears to be a reimported in-progress project. The design assets are present and the main workflow is intact:

- `SKILL.md`: repo-level operating guide for DSM design work
- `references/designs/dt/dsm_quick_design.m`: editable quick-start template
- `references/designs/dt/design_3rd_order_ciff_10mhz.m`: 3rd-order lowpass example with amplitude sweep
- `references/designs/ct/design_3rd_order_ct_dsm_10mhz.m`: 3rd-order continuous-time mapping example using `realizeNTF_ct`
- `references/designs/dt/design_4th_order_ciff.m`: full 4th-order CIFF example
- `references/designs/dt/design_4th_order_ciff_with_dac_mismatch.m`: mismatch and DWA study
- `references/simulink/builders/build_dsm_simulink_model.m`: generated Simulink builder from DSM coefficients or ABCD data
- `references/simulink/builders/build_ct_dsm_simulink_model.m`: generic CTDSM Simulink builder with CT integrators, `fs` scaling, and sampled quantizer input
- `references/analysis/jitter/calculate_ct_dac_jitter_sjnr.m`: standalone `J` and `SJNR` calculator for CT DAC clock jitter
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_simulink_model.m`: CT DAC clock-jitter sweep for the compensated 3rd-order CTDSM
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_theory_compare.m`: compare Simulink jitter modes against the equation-based SJNR estimate
- `references/simulink/runs/run_3rd_order_ct_eld_simulink_model.m`: Ts/2 excess-loop-delay comparison for uncompensated and compensated CTDSM cases
- `references/simulink/runs/run_3rd_order_ct_simulink_model.m`: CT Simulink execution and plotting flow
- `references/simulink/runs/run_3rd_order_simulink_model.m`: Simulink execution and plotting flow
- `references/simulink/tests/test_build_3rd_order_ct_simulink_model.m`: smoke test for the generic CT Simulink builder
- `references/third_party/dstoolbox/`: bundled Delta Sigma Toolbox sources

There are also saved outputs in `references/*.mat`, `references/*.png`, and `references/*.slx`.

## Environment Notes

The local toolchain has now been repaired:

- Apple Command Line Tools are installed and active.
- `git` and `python3` work from the shell again.
- MATLAB is installed at `/Applications/MATLAB_R2025b.app`.
- A user-level `matlab` launcher is available from normal `zsh` sessions via `~/bin/matlab`.

If a tool does not resolve inside a long-lived app shell, open a fresh terminal or start a new login shell so the updated `PATH` is reloaded.

## Verified Baseline

The 3rd-order reference flow has been rerun successfully on April 12, 2026:

- Direct MATLAB script: `references/designs/dt/design_3rd_order_ciff_10mhz.m`
- Continuous-time mapped script: `references/designs/ct/design_3rd_order_ct_dsm_10mhz.m`
- Continuous-time Simulink build test: `references/simulink/tests/test_build_3rd_order_ct_simulink_model.m`
- Continuous-time Simulink run: `references/simulink/runs/run_3rd_order_ct_simulink_model.m`
- Continuous-time Simulink ELD comparison: `references/simulink/runs/run_3rd_order_ct_eld_simulink_model.m`
- Generated Simulink flow: `references/simulink/runs/run_3rd_order_simulink_model.m`

Saved results:

- `references/results/dt/dsm_3rd_order_4bit_osr32_results.mat`
- `references/results/dt/dsm_3rd_order_4bit_osr32_plots.png`
- `references/results/ct/base/dsm_3rd_order_ct_10mhz_results.mat`
- `references/results/ct/base/dsm_3rd_order_ct_10mhz_plots.png`
- `references/results/simulink/ct/dsm_3rd_order_ct_10mhz_topology_model.slx`
- `references/results/simulink/ct/dsm_3rd_order_ct_10mhz_simulink_results.mat`
- `references/results/simulink/ct/dsm_3rd_order_ct_10mhz_simulink_plots.png`
- `references/results/ct/eld/dsm_3rd_order_ct_eld_simulink_results.mat`
- `references/results/ct/eld/dsm_3rd_order_ct_eld_simulink_plots.png`
- `references/results/ct/jitter/dsm_3rd_order_ct_dac_jitter_simulink_results.mat`
- `references/results/ct/jitter/dsm_3rd_order_ct_dac_jitter_simulink_plots.png`
- `references/results/ct/jitter/dsm_3rd_order_ct_dac_jitter_theory_compare.mat`
- `references/results/ct/jitter/dsm_3rd_order_ct_dac_jitter_theory_compare.png`
- `references/results/simulink/dt/dsm_3rd_order_ciff_10mhz_simulink_results.mat`
- `references/results/simulink/dt/dsm_3rd_order_ciff_10mhz_simulink_plots.png`

Measured baseline:

- Direct MATLAB: `SNR = 109.27 dB`, `ENOB = 17.86 bits`
- Direct MATLAB sweep peak: `112.44 dB` at `0.85 V`
- Continuous-time mapped DT-equivalent: `SNR = 83.06 dB`, `ENOB = 13.51 bits`
- Continuous-time mapped sweep peak: `89.41 dB` at `0.90 V`
- Continuous-time Simulink: `SNR = 106.97 dB`, `ENOB = 17.48 bits`
- Continuous-time Simulink with `Ts/2` ELD, uncompensated: `SNR = 27.74 dB`, `ENOB = 4.32 bits`
- Continuous-time Simulink with `Ts/2` ELD, compensated: `SNR = 107.22 dB`, `ENOB = 17.52 bits`
- Continuous-time Simulink with `Ts/2` ELD and `500 ps` DAC edge jitter: `SNR = 66.02 dB`
- Continuous-time Simulink with `Ts/2` ELD and `500 ps` equivalent DAC jitter: `SNR = 66.73 dB`
- Equation-style CT DAC jitter estimate for the same `500 ps` case: `J = 2.98e-8`, `SJNR = 66.22 dB`
- Simulink: `SNR = 109.21 dB`, `ENOB = 17.85 bits`

The direct MATLAB and Simulink paths agree within about `0.06 dB`, which is a strong indication that the current 3rd-order design path is intact.

The current CT path uses the toolbox `FF` continuous-time realization so the generated Simulink topology stays CIFF-like, keeps an explicit sampler ahead of the quantizer, and scales each CT integrator drive by `fs` for correct normalization.

The CT builder now also supports delayed DAC pulses with width up to one sample period, which allows direct modeling of excess loop delay. For the verified half-cycle ELD case, the uncompensated model uses the original CT coefficients with delayed DAC timing `tdac = [0.5 1.5]`, while the compensated model uses `realizeNTF_ct(ntf, 'FF', [0.5 1.5])` so the compensation path is synthesized automatically by the toolbox.

The CT builder also supports two DAC clock-jitter modes for full-width NRZ-like delayed DAC pulses:

- `jitter_mode='edge'`: perturbs the DAC update edge in continuous time
- `jitter_mode='equivalent'`: injects the textbook-style sampled jitter error term derived from `v[n] - v[n-1]`

For the current 3rd-order compensated CTDSM, both jitter modes produce similar SNR at `500 ps` RMS, and the equation-style estimate is now close as well once the NTF integral is converted from the numerically evaluated `f` domain to the `w` domain used by Eq. 9.17 via a `2*pi` factor. When doing coherent tone measurements, remember MATLAB FFT arrays are 1-based: a tone at DFT bin `k` appears at array index `k+1`.

## Good Resume Points

If you want to continue the project, these are the most likely starting points:

1. Run `references/designs/dt/dsm_quick_design.m` for a fresh parameterized DSM design.
2. Run `references/designs/dt/design_3rd_order_ciff_10mhz.m` to validate the current 3rd-order reference flow.
3. Run `references/designs/ct/design_3rd_order_ct_dsm_10mhz.m` to validate the continuous-time mapping flow.
4. Run `references/simulink/tests/test_build_3rd_order_ct_simulink_model.m` to build and smoke-test the CT Simulink model.
5. Run `references/simulink/runs/run_3rd_order_ct_simulink_model.m` to simulate and plot the CT Simulink model.
6. Run `references/simulink/runs/run_3rd_order_ct_eld_simulink_model.m` to compare uncompensated and compensated `Ts/2` excess-loop-delay behavior.
7. Run `references/simulink/runs/run_3rd_order_ct_dac_jitter_simulink_model.m` to sweep DAC clock jitter in the CT Simulink model.
8. Run `references/analysis/jitter/calculate_ct_dac_jitter_sjnr.m` to reproduce the current equation-style `J` and `SJNR` calculation.
9. Run `references/simulink/runs/run_3rd_order_ct_dac_jitter_theory_compare.m` to compare the Simulink jitter modes against the equation-based estimate.
10. Run `references/simulink/runs/run_3rd_order_simulink_model.m` to verify the generated discrete-time Simulink topology path.
11. Run `references/designs/dt/design_4th_order_ciff_with_dac_mismatch.m` if the next focus is DAC mismatch or DWA behavior.

## Recovery Work Done

During recovery, stale hardcoded paths from the original workspace were identified and removed from the repo files so the project is portable again after reimport.
