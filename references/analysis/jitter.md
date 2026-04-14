# CT DAC Jitter

Use this guide when the task involves DAC clock jitter in a CT DSM.

## Modes

- `jitter_mode='edge'`: perturb the DAC update edge directly in Simulink for waveform-level jitter studies.
- `jitter_mode='equivalent'`: inject a sampled error term proportional to `v[n] - v[n-1]` for textbook-style comparisons.

## Modeling Notes

- For delayed NRZ-like CT feedback pulses, start with `jitter_mode='edge'`.
- Use `jitter_mode='equivalent'` when comparing against equation-based SJNR estimates.
- Ensure the DAC pulse timing in simulation matches the `tdac` used during CT realization.

## Primary References

- `references/analysis/jitter/calculate_ct_dac_jitter_sjnr.m`: standalone equation-style `J` and `SJNR` calculation.
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_simulink_model.m`: simulate DAC clock-jitter sweeps.
- `references/simulink/runs/run_3rd_order_ct_dac_jitter_theory_compare.m`: compare edge-jitter, equivalent-jitter, and theory side by side.
