# Codex_DSM

MATLAB reference workspace for delta-sigma modulator design using Richard Schreier's Delta Sigma Toolbox, plus supporting examples for flash ADC quantization, thermometer DAC feedback, DAC mismatch, DWA, and generated Simulink models.

## Current State

This repo appears to be a reimported in-progress project. The design assets are present and the main workflow is intact:

- `SKILL.md`: repo-level operating guide for DSM design work
- `references/dsm_quick_design.m`: editable quick-start template
- `references/design_3rd_order_ciff_10mhz.m`: 3rd-order lowpass example with amplitude sweep
- `references/design_3rd_order_ct_dsm_10mhz.m`: 3rd-order continuous-time mapping example using `realizeNTF_ct`
- `references/design_4th_order_ciff.m`: full 4th-order CIFF example
- `references/design_4th_order_ciff_with_dac_mismatch.m`: mismatch and DWA study
- `references/build_dsm_simulink_model.m`: generated Simulink builder from DSM coefficients or ABCD data
- `references/build_ct_dsm_simulink_model.m`: generic CTDSM Simulink builder with CT integrators, `fs` scaling, and sampled quantizer input
- `references/run_3rd_order_ct_eld_simulink_model.m`: Ts/2 excess-loop-delay comparison for uncompensated and compensated CTDSM cases
- `references/run_3rd_order_ct_simulink_model.m`: CT Simulink execution and plotting flow
- `references/run_3rd_order_simulink_model.m`: Simulink execution and plotting flow
- `references/test_build_3rd_order_ct_simulink_model.m`: smoke test for the generic CT Simulink builder
- `references/dstoolbox/`: bundled Delta Sigma Toolbox sources

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

- Direct MATLAB script: `references/design_3rd_order_ciff_10mhz.m`
- Continuous-time mapped script: `references/design_3rd_order_ct_dsm_10mhz.m`
- Continuous-time Simulink build test: `references/test_build_3rd_order_ct_simulink_model.m`
- Continuous-time Simulink run: `references/run_3rd_order_ct_simulink_model.m`
- Continuous-time Simulink ELD comparison: `references/run_3rd_order_ct_eld_simulink_model.m`
- Generated Simulink flow: `references/run_3rd_order_simulink_model.m`

Saved results:

- `references/dsm_3rd_order_4bit_osr32_results.mat`
- `references/dsm_3rd_order_4bit_osr32_plots.png`
- `references/dsm_3rd_order_ct_10mhz_results.mat`
- `references/dsm_3rd_order_ct_10mhz_plots.png`
- `references/dsm_3rd_order_ct_10mhz_topology_model.slx`
- `references/dsm_3rd_order_ct_10mhz_simulink_results.mat`
- `references/dsm_3rd_order_ct_10mhz_simulink_plots.png`
- `references/dsm_3rd_order_ct_eld_simulink_results.mat`
- `references/dsm_3rd_order_ct_eld_simulink_plots.png`
- `references/dsm_3rd_order_ciff_10mhz_simulink_results.mat`
- `references/dsm_3rd_order_ciff_10mhz_simulink_plots.png`

Measured baseline:

- Direct MATLAB: `SNR = 109.27 dB`, `ENOB = 17.86 bits`
- Direct MATLAB sweep peak: `112.44 dB` at `0.85 V`
- Continuous-time mapped DT-equivalent: `SNR = 83.06 dB`, `ENOB = 13.51 bits`
- Continuous-time mapped sweep peak: `89.41 dB` at `0.90 V`
- Continuous-time Simulink: `SNR = 106.97 dB`, `ENOB = 17.48 bits`
- Continuous-time Simulink with `Ts/2` ELD, uncompensated: `SNR = 27.74 dB`, `ENOB = 4.32 bits`
- Continuous-time Simulink with `Ts/2` ELD, compensated: `SNR = 107.22 dB`, `ENOB = 17.52 bits`
- Simulink: `SNR = 109.21 dB`, `ENOB = 17.85 bits`

The direct MATLAB and Simulink paths agree within about `0.06 dB`, which is a strong indication that the current 3rd-order design path is intact.

The current CT path uses the toolbox `FF` continuous-time realization so the generated Simulink topology stays CIFF-like, keeps an explicit sampler ahead of the quantizer, and scales each CT integrator drive by `fs` for correct normalization.

The CT builder now also supports delayed DAC pulses with width up to one sample period, which allows direct modeling of excess loop delay. For the verified half-cycle ELD case, the uncompensated model uses the original CT coefficients with delayed DAC timing `tdac = [0.5 1.5]`, while the compensated model uses `realizeNTF_ct(ntf, 'FF', [0.5 1.5])` so the compensation path is synthesized automatically by the toolbox.

## Good Resume Points

If you want to continue the project, these are the most likely starting points:

1. Run `references/dsm_quick_design.m` for a fresh parameterized DSM design.
2. Run `references/design_3rd_order_ciff_10mhz.m` to validate the current 3rd-order reference flow.
3. Run `references/design_3rd_order_ct_dsm_10mhz.m` to validate the continuous-time mapping flow.
4. Run `references/test_build_3rd_order_ct_simulink_model.m` to build and smoke-test the CT Simulink model.
5. Run `references/run_3rd_order_ct_simulink_model.m` to simulate and plot the CT Simulink model.
6. Run `references/run_3rd_order_ct_eld_simulink_model.m` to compare uncompensated and compensated `Ts/2` excess-loop-delay behavior.
7. Run `references/run_3rd_order_simulink_model.m` to verify the generated discrete-time Simulink topology path.
8. Run `references/design_4th_order_ciff_with_dac_mismatch.m` if the next focus is DAC mismatch or DWA behavior.

## Recovery Work Done

During recovery, stale hardcoded paths from the original workspace were identified and removed from the repo files so the project is portable again after reimport.
