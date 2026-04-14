# MASH Design

Use this guide for discrete-time cascaded residue-shaping MASH work.

## Core Sequence

1. Design or provide one NTF per stage.
2. Realize each stage with `realizeNTF(..., 'CIFF')` or the requested form.
3. Drive stage 1 from the external input and each later stage from the previous stage residue.
4. Form the residue as `y_k - v_k`.
5. For textbook cascaded residue-shaping MASH cancellation, combine the stage outputs digitally with cumulative previous-stage NTF products:
   `v_mash = v_1 + NTF_1(z)v_2 + NTF_1(z)NTF_2(z)v_3 + ...`
6. Prefer explicit stage-by-stage simulation when debugging residue transfer or cancellation.

## Defaults

| Parameter | Typical default | Notes |
| --- | --- | --- |
| `MASH structure` | `2-1` | Good first cascaded DT reference example |

## Guidance

- `MASH 2-1`: good introductory DT cascaded-noise-shaping case when you want a 3rd-order overall response from stable lower-order stages.
- Use explicit stage outputs and residues when debugging cancellation or unexpected tone leakage.

## Primary References

- `references/designs/dt/design_mash_2_1_ciff_10mhz.m`: best starting point for the current DT 2-1 MASH reference flow.
- `references/simulink/builders/build_mash_dsm_simulink_model.m`: generic DT MASH Simulink builder from stage descriptions and digital cancellation filters.
- `references/simulink/builders/build_mash_21_dsm_simulink_model.m`: thin compatibility wrapper for the current 2-1 MASH example.
- `references/simulink/runs/run_mash_2_1_simulink_model.m`: execute and plot the current MASH Simulink model.
