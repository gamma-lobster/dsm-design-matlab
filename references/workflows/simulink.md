# Simulink Workflows

Use this guide when the user wants the design assembled or exercised in Simulink.

## DT And MASH Builders

- `references/simulink/builders/build_dsm_simulink_model.m`: assemble a Simulink model from a DT DSM design.
- `references/simulink/builders/build_mash_dsm_simulink_model.m`: assemble a generic DT MASH Simulink model from stage descriptions and digital cancellation filters.
- `references/simulink/builders/build_mash_21_dsm_simulink_model.m`: compatibility wrapper for the current 2-1 MASH example.

## CT Builders

- `references/simulink/builders/build_ct_dsm_simulink_model.m`: assemble a generic CTDSM Simulink model with CT integrators and a sampled quantizer input.

## Runs

- `references/simulink/runs/run_mash_2_1_simulink_model.m`: run and plot the current MASH Simulink model.
- `references/simulink/runs/run_3rd_order_ct_simulink_model.m`: run and plot the CT Simulink model.
- `references/simulink/runs/run_3rd_order_ct_eld_simulink_model.m`: compare uncompensated and compensated ELD behavior.
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_simulink_model.m`: run a DAC clock-jitter sweep in the CT Simulink model.
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_theory_compare.m`: compare simulated and equation-style jitter results side by side.

## Simulink Rules

- For CT models, scale integrator drives by `fs` so normalized toolbox coefficients map correctly into seconds.
- Place an explicit sampler before the ADC or quantizer in CT models.
- Ensure Simulink DAC timing matches the `tdac` used during `realizeNTF_ct(...)`.
