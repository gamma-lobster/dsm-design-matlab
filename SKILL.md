---
name: dsm-design-matlab
description: Design and analyze delta-sigma modulators in MATLAB using Richard Schreier's Delta Sigma Toolbox. Use when Codex needs to synthesize an NTF, realize DT or CT DSM topologies, build ABCD or ABCDc matrices, assemble MASH structures, simulate SNR or ENOB, troubleshoot stability, or model multi-bit flash ADC and thermometer DAC behavior for DSM work.
---

# DSM Design In MATLAB

Use the native MATLAB Delta Sigma Toolbox workflow and keep the top-level skill lean. Open the scenario-specific Markdown files under `references/` only when they match the user request.

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

For discrete-time DSM work, use:

```matlab
synthesizeNTF -> realizeNTF -> stuffABCD -> simulate/analyze
```

For continuous-time DSM work, use:

```matlab
synthesizeNTF -> realizeNTF_ct -> mapCtoD -> simulate/analyze
```

Keep designs normalized unless the user gives physical circuit scaling requirements.

## Use The Bundled References

Open only the reference files needed for the current task. Keep reference loading shallow: `SKILL.md` should be the dispatcher, and the scenario docs below should be the first place to branch.

### Workflow Guides

- `references/workflows/dt_design.md`: standard discrete-time DSM design flow, defaults, topology guidance, and starter files.
- `references/workflows/ct_design.md`: continuous-time DSM design, `realizeNTF_ct`, `mapCtoD`, ELD, and jitter guidance.
- `references/workflows/mash.md`: cascaded DT MASH structure, residue flow, and cancellation guidance.
- `references/workflows/simulink.md`: which Simulink builders and runners to use for DT, MASH, and CT tasks.

### Analysis Guides

- `references/analysis/snr_enob.md`: FFT-based SNR, SNDR, and ENOB rules.
- `references/analysis/stability_debug.md`: stability and debugging checklist.
- `references/analysis/jitter.md`: CT DAC jitter modes and related references.

### Component Guides

- `references/components/multi_bit_feedback.md`: flash ADC, thermometer DAC, mismatch, and DWA guidance.

If the user asks about toolbox behavior, inspect the relevant source in `references/third_party/dstoolbox/` rather than guessing.

## Response Pattern

When helping with a design, prefer this output structure:

1. Restate the target spec in MATLAB terms.
2. Show the chosen design parameters and why they are reasonable.
3. Provide or edit MATLAB code that follows the DT or CT toolbox workflow.
4. Summarize the expected stability and performance limits.
5. Call out any assumptions, especially around OSR, topology, stage partitioning for MASH, quantizer bits, DAC timing, ELD compensation, and mismatch.

Keep explanations practical and tied to the provided MATLAB files.
