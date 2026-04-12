%% calculate_ct_dac_jitter_sjnr.m
% Compute CTDSM DAC-jitter noise power J and jitter-limited SJNR using
% the f-domain integral
%
%   int_0^(1/2) |(1 - z^-1)NTF(z)|^2 df
%
% where:
%   - f is normalized discrete-time frequency in cycles/sample
%   - z = exp(j*2*pi*f)
%
% This script evaluates the integral over the full normalized band
% 0 <= f <= 1/2, then applies the Eq. 9.17-style scaling used in this
% project:
%
%   J = (sigma_t^2/T^2) * (1/OSR) * (1/(3*pi)) * int_term * (Delta/2)^2
%
%   SJNR = 10*log10(P_sig/J)
%
% Assumptions used here:
%   - single DAC clock jitter RMS value sigma_t
%   - quantizer/DAC step size Delta = 2*V_fs/(n_levels - 1)
%   - input signal power P_sig = A_in^2/2 for a sinusoidal input
%
% The defaults match the current 3rd-order CTDSM reference design.

clearvars; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(script_dir, 'dstoolbox'));

order = 3;
OSR = 32;
opt = 1;
H_inf = 4.0;
f0 = 0;
n_grid = 200000;

fs = 10e6;
T = 1 / fs;
sigma_t = 500e-12;

n_bits = 4;
n_levels = 2^n_bits;
V_fs = 1.0;
A_in = 0.5;

ntf = synthesizeNTF(order, OSR, opt, H_inf, f0);

% Evaluate the DT NTF on the unit circle using normalized frequency f.
f = linspace(0, 0.5, n_grid);
z = exp(1j * 2*pi * f);
H = evalTF(ntf, z);

% Full-band f-domain integral of the DAC-jitter weighting term.
integrand = abs((1 - z.^(-1)) .* H).^2;
int_term = trapz(f, integrand);

% Actual multibit DAC step size for the quantizer used in this repo.
Delta = 2 * V_fs / (n_levels - 1);

% Input sine-wave power for amplitude A_in.
P_sig = A_in^2 / 2;

% Eq. 9.17-style jitter noise power and corresponding jitter-limited SJNR.
J = (sigma_t^2 / T^2) * (1 / OSR) * (1 / (3*pi)) * int_term * (Delta / 2)^2;
SJNR = 10 * log10(P_sig / J);

fprintf('CTDSM DAC jitter calculation\n');
fprintf('  order      = %d\n', order);
fprintf('  OSR        = %d\n', OSR);
fprintf('  H_inf      = %.6g\n', H_inf);
fprintf('  fs         = %.6g Hz\n', fs);
fprintf('  sigma_t    = %.6g s\n', sigma_t);
fprintf('  n_levels   = %d\n', n_levels);
fprintf('  V_fs       = %.6g\n', V_fs);
fprintf('  A_in       = %.6g\n', A_in);
fprintf('  grid pts   = %d\n', n_grid);
fprintf('  band       = [0, 0.5] cycles/sample\n');
fprintf('  int_term   = %.12e\n', int_term);
fprintf('  Delta      = %.12e\n', Delta);
fprintf('  P_sig      = %.12e\n', P_sig);
fprintf('  J          = %.12e\n', J);
fprintf('  SJNR       = %.12f dB\n', SJNR);

results = struct();
results.order = order;
results.OSR = OSR;
results.opt = opt;
results.H_inf = H_inf;
results.f0 = f0;
results.n_grid = n_grid;
results.fs = fs;
results.T = T;
results.sigma_t = sigma_t;
results.n_bits = n_bits;
results.n_levels = n_levels;
results.V_fs = V_fs;
results.A_in = A_in;
results.f = f;
results.z = z;
results.H = H;
results.integrand = integrand;
results.int_term = int_term;
results.Delta = Delta;
results.P_sig = P_sig;
results.J = J;
results.SJNR = SJNR;

save(fullfile(script_dir, 'ct_dac_jitter_sjnr_results.mat'), 'results');
