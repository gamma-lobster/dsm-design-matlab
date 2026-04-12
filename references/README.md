# DSM Design Examples

This folder contains MATLAB reference designs and the Delta Sigma Toolbox.

## Folder Structure

```
references/
├── dstoolbox/              # Delta Sigma Toolbox (117 functions)
│   ├── synthesizeNTF.m
│   ├── realizeNTF.m
│   ├── stuffABCD.m
│   ├── simulateDSM.m
│   └── ... (113 more)
│
├── design_3rd_order_ciff_10mhz.m # 3rd-order, 4-bit, OSR=32, fs=10MHz example
├── design_4th_order_ciff.m    # Full 4th-order CIFF example
├── dsm_4th_order_simple.m     # Simplified headless version
├── dsm_quick_design.m         # Quick-start template
│
├── dsm_4th_order_results.mat  # Simulation results
├── dsm_plots.html           # HTML plots
└── create_plots.py          # Python plot generator
```

## Quick Start

All scripts automatically add the toolbox to MATLAB's path:

```matlab
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));
```

Just run any `.m` file directly - no manual path setup needed!

## Example: 4th-Order CIFF

```matlab
% From MATLAB
cd /path/to/references
run('design_4th_order_ciff.m')
```

## Results

| Metric | Value |
|--------|-------|
| SNR | 103.72 dB |
| ENOB | 16.94 bits |
| Status | Stable |

## Verified 3rd-Order Baseline

The 3rd-order 10 MHz example and its generated Simulink model were rerun successfully on April 12, 2026.

| Flow | SNR | ENOB | Notes |
|------|-----|------|-------|
| `design_3rd_order_ciff_10mhz.m` | 109.27 dB | 17.86 bits | Peak sweep SNR 112.44 dB at 0.85 V |
| `run_3rd_order_simulink_model.m` | 109.21 dB | 17.85 bits | Simulink output matches MATLAB within about 0.06 dB |

Fresh outputs from that verification:

- `dsm_3rd_order_4bit_osr32_results.mat`
- `dsm_3rd_order_4bit_osr32_plots.png`
- `dsm_3rd_order_ciff_10mhz_simulink_results.mat`
- `dsm_3rd_order_ciff_10mhz_simulink_plots.png`

## Notes

- `dstoolbox/` contains Richard Schreier's Delta Sigma Toolbox
- All scripts use **native MATLAB only** - no Python dependencies
- Based on "Understanding Delta-Sigma Data Converters" (Appendix B)
