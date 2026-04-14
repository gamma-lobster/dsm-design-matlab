# DSM Design Examples

This folder is organized by role so the MATLAB references are easier to
navigate and extend.

## Structure

```text
references/
├── common/                 shared MATLAB helpers
├── designs/                direct MATLAB design flows
│   ├── dt/
│   └── ct/
├── simulink/               generated-model builders, runs, and tests
│   ├── builders/
│   ├── runs/
│   └── tests/
├── analysis/               standalone analysis and debug scripts
│   ├── jitter/
│   └── debug/
├── components/             flash ADC, thermometer DAC, mismatch, DWA helpers
├── results/                saved plots, MAT files, and generated models
│   ├── dt/
│   ├── ct/
│   └── simulink/
├── tools/                  optional Python helpers and extracted assets
└── third_party/dstoolbox/  Richard Schreier's Delta Sigma Toolbox
```

## Quick Start

The MATLAB entry points now use `common/setup_references_paths.m` to find
the reorganized toolbox, component, and builder folders automatically.

Useful starting points:

- `designs/dt/dsm_quick_design.m`
- `designs/dt/design_3rd_order_ciff_10mhz.m`
- `designs/ct/design_3rd_order_ct_dsm_10mhz.m`
- `simulink/runs/run_3rd_order_ct_simulink_model.m`
- `analysis/jitter/calculate_ct_dac_jitter_sjnr.m`

The skill guidance is also split by scenario so Codex can load only what it needs:

- `../SKILL.md`: top-level dispatcher
- `workflows/dt_design.md`: discrete-time design flow
- `workflows/ct_design.md`: continuous-time design, ELD, and jitter flow
- `workflows/mash.md`: cascaded MASH flow
- `workflows/simulink.md`: Simulink builders and runs
- `analysis/snr_enob.md`: FFT-based measurement rules
- `analysis/stability_debug.md`: debugging checklist
- `components/multi_bit_feedback.md`: flash ADC, DAC mismatch, and DWA guidance

## Results Layout

Generated outputs are grouped by topic instead of being mixed with source:

- `results/dt/`: direct MATLAB DT design outputs
- `results/ct/base/`: direct MATLAB CT mapping outputs
- `results/ct/eld/`: excess-loop-delay runs
- `results/ct/jitter/`: DAC clock-jitter runs and analysis
- `results/simulink/dt/`: DT Simulink models and outputs
- `results/simulink/ct/`: CT Simulink models and outputs

## Notes

- `third_party/dstoolbox/` contains Richard Schreier's Delta Sigma Toolbox.
- The MATLAB scripts are still intended to run directly without manual path setup.
- `tools/` contains optional helpers; the main design flows remain native MATLAB.
