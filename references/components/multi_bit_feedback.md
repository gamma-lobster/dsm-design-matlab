# Multi-Bit Feedback Components

Use this guide when the DSM task includes a flash ADC, thermometer DAC reconstruction, DAC mismatch, or DWA.

## Supported Tasks

- model a flash ADC that outputs thermometer code
- reconstruct the analog level with a thermometer DAC
- inject unit-cell mismatch into a multi-bit DAC
- compare ideal, mismatched, and DWA-shaped feedback DAC behavior

Use the bundled component functions before writing new helpers. Keep thermometer outputs contiguous unless the task explicitly requires rotation or scrambling. For DWA, rotate cell usage for DAC feedback while preserving the ADC's original code meaning.

## Primary References

- `references/components/flash_adc_quantizer.m`: flash ADC with thermometer-code output.
- `references/components/thermometer_dac.m`: ideal thermometer DAC reconstruction.
- `references/components/thermometer_dac_mismatch.m`: static unit-cell mismatch model.
- `references/components/thermometer_dac_dwa.m`: dynamic weighted averaging for mismatch shaping.
- `references/designs/dt/design_4th_order_ciff_with_dac_mismatch.m`: end-to-end mismatch and DWA example.
