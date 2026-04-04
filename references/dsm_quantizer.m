function v = dsm_quantizer(y, V_fs, n_levels)
%DSM_QUANTIZER Ideal symmetric mid-tread/mid-rise quantizer for Simulink.

lsb = 2 * V_fs / (n_levels - 1);
levels = linspace(-V_fs, V_fs, n_levels);
y_clip = min(max(y, -V_fs), V_fs);
idx = round((y_clip + V_fs) / lsb) + 1;
idx = min(max(idx, 1), n_levels);
v = levels(idx);
end
