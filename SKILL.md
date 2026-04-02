---
name: dsm-design-matlab
description: "Delta-Sigma modulator design using native MATLAB Delta Sigma Toolbox. Pure MATLAB workflow: synthesizeNTF, realizeNTF, stuffABCD, and simulation. No Python dependency. Based on Richard Schreier's methodology from 'Understanding Delta-Sigma Data Converters'."
homepage: https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox
metadata:
  {
    "openclaw":
      {
        "emoji": "🔺",
        "requires": { "bins": ["matlab"] },
        "optional": {},
        "install": [],
      },
  }
---

# Delta-Sigma Modulator Design Skill (MATLAB Native)

This skill guides you through designing Delta-Sigma modulators using the native MATLAB Delta Sigma Toolbox by Richard Schreier, as documented in Appendix B of "Understanding Delta-Sigma Data Converters" (2nd Edition) by Pavan, Schreier, and Temes.

## MATLAB Native Workflow

Unlike the Octave version which requires Python for NTF synthesis, MATLAB can run the full workflow natively:

```
MATLAB: synthesizeNTF → realizeNTF → stuffABCD → simulateDSM
```

## Prerequisites

### Delta Sigma Toolbox Installation

Download from MATLAB Central: https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox

**Recommended: Place toolbox in a subfolder (e.g., `dstoolbox/`)**

```matlab
% Add toolbox to path (relative to script location)
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));
```

Or for a global installation:
```matlab
% Add toolbox to path
addpath('/usr/local/MATLAB/toolbox/deltasig');
savepath;  % Optional: save for future sessions
```

### Verify Installation

```matlab
% Test basic functions
which synthesizeNTF
which realizeNTF
which stuffABCD
which simulateDSM
```

## Quick Start: Native MATLAB Workflow

### Step 1: Synthesize NTF

```matlab
%% NTF Synthesis
order = 3;      % Modulator order
OSR = 64;       % Oversampling ratio
H_inf = 1.5;    % Max out-of-band gain (Lee's rule: < 2.0)
opt = 1;        % Optimize zeros (0 = no opt, 1 = opt)
f0 = 0;         % Center frequency (0 = LP, 0.25 = BP at fs/4)

% Synthesize NTF
[ntf, stf] = synthesizeNTF(order, OSR, opt, H_inf, f0);

% Display
zpk(ntf)
```

### Step 2: Realize Coefficients

```matlab
%% Realize NTF to coefficients
form = 'CIFF';  % Cascade of Integrators, Feedback Form
                % Options: 'CIFF', 'CIFB', 'CRFB', 'CRFF', 'CRFBD'

[a, g, b, c] = realizeNTF(ntf, form);

fprintf('Coefficients:\n');
fprintf('  a = ['); fprintf('%.6f ', a); fprintf('] (feedback)\n');
fprintf('  g = ['); fprintf('%.6f ', g); fprintf('] (resonator)\n');
fprintf('  b = ['); fprintf('%.6f ', b); fprintf('] (input)\n');
fprintf('  c = ['); fprintf('%.6f ', c); fprintf('] (inter-stage)\n');
```

### Step 3: Build ABCD Matrix

```matlab
%% Build ABCD using stuffABCD
ABCD = stuffABCD(a, g, b, c, form);

% Verify
fprintf('\nABCD Matrix:\n');
disp(ABCD);

% Extract for simulation
[n_states, ~] = size(ABCD);
n_states = n_states - 1;  % Last row is output
A_mat = ABCD(1:n_states, 1:n_states);
B_mat = ABCD(1:n_states, n_states+1:end);
C_mat = ABCD(n_states+1, 1:n_states);
D_mat = ABCD(n_states+1, n_states+1:end);
```

### Step 4: Simulate

```matlab
%% Time-domain simulation
fs = 1e6;           % Sampling frequency (Hz)
N = 8192;           % Number of samples
fB = fs / (2*OSR);  % Signal bandwidth

% Input signal
f_bin = round(sqrt(1/7) * N / (2*OSR));  % Not at harmonic
A_in = 0.5;         % Input amplitude
u = A_in * sin(2*pi*f_bin*(0:N-1)/N);

% Quantizer
n_bits = 5;
n_levels = 2^n_bits;
V_fs = 1.0;
LSB = 2*V_fs / (n_levels - 1);
q_levels = linspace(-V_fs, V_fs, n_levels);

% Run simulation
x = zeros(n_states, 1);
y = zeros(1, N);
v = zeros(1, N);

for i = 1:N
    if i == 1, v_prev = 0; else, v_prev = v(i-1); end
    
    % Compute output
    y(i) = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;
    
    % Quantize
    y_clip = max(-V_fs, min(V_fs, y(i)));
    idx = round((y_clip - (-V_fs)) / LSB) + 1;
    idx = max(1, min(n_levels, idx));
    v(i) = q_levels(idx);
    
    % Update state
    x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v(i);
end
```

### Step 5: Calculate SNR

```matlab
%% Calculate SNR from output spectrum
w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));  % Hann window
V_out = fft(v .* w) / (N/4);
V_out_mag = abs(V_out);

% Find signal
[~, sig_idx] = max(V_out_mag(2:N/2));
sig_bin = sig_idx + 1;
fB_bins = ceil(N / (2*OSR));

% Signal bins
sig_bins = sig_bin-1:sig_bin+1;
sig_bins = sig_bins(sig_bins >= 2 & sig_bins <= fB_bins);
signal_power = sum(V_out_mag(sig_bins).^2);

% Noise bins
noise_bins = setdiff(2:fB_bins, sig_bins);
noise_power = sum(V_out_mag(noise_bins).^2);

% Calculate SNR
SNR = 10*log10(signal_power / noise_power);
ENOB = (SNR - 1.76) / 6.02;

fprintf('SNR: %.2f dB\n', SNR);
fprintf('ENOB: %.2f bits\n', ENOB);
```

## Core Concepts

### Key Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `order` | Modulator order | 1-5 for stable designs |
| `OSR` | Oversampling Ratio | 32, 64, 128 |
| `H_inf` | Max out-of-band NTF gain | 1.5-2.0 (Lee's rule: < 2) |
| `f0` | Center frequency | 0 for LP, 0.25 for BP at fs/4 |
| `opt` | Zero optimization | 0=none, 1=optimized |
| `form` | Topology | 'CIFF', 'CIFB', 'CRFB', 'CRFF' |
| `n_bits` | Quantizer bits | 1-8 |
| `n_levels` | Quantizer levels | 2^n_bits |
| `V_fs` | Full scale voltage | +/- 1 (normalized) |

### Topology Options

| Topology | Form | Characteristics |
|----------|------|-----------------|
| CIFF | Cascade of Integrators, Feedback | Unity STF, good for DAC feedback |
| CIFB | Cascade of Integrators, FeedBack | Feedforward summation |
| CRFB | Cascade of Resonators, FeedBack | For bandpass designs |
| CRFF | Cascade of Resonators, FeedForward | Feedforward resonator outputs |
| CRFBD | Cascade of Resonators, FeedBack with feedforward | Complex pole NTF |

### Design Flow (MATLAB Native)

```
┌─────────────────────────────────────────────────────────────┐
│  MATLAB: Complete Native Workflow                           │
│  ────────────────────────────────                           │
│                                                             │
│  1. Call synthesizeNTF(order, OSR, opt, H_inf, f0)        │
│     Returns: ntf (zpk object), stf                          │
│                                                             │
│  2. Call realizeNTF(ntf, 'CIFF')                          │
│     Returns: a, g, b, c coefficients                      │
│                                                             │
│  3. Call stuffABCD(a, g, b, c, 'CIFF')                      │
│     Returns: ABCD state-space matrix                        │
│                                                             │
│  4. Run time-domain simulation                              │
│     State update: x = A*x + B*[u; v]                        │
│     Output: y = C*x + D*[u; v]                              │
│                                                             │
│  5. Calculate SNR from output spectrum                      │
│     Find signal in OUTPUT spectrum (not noise spectrum)   │
│     ENOB = (SNR - 1.76) / 6.02                              │
│                                                             │
│  6. Plot NTF, STF, spectra                                  │
│     Use toolbox plotNTF() or custom plots                   │
└─────────────────────────────────────────────────────────────┘
```

## Useful MATLAB Functions from Delta Sigma Toolbox

### Core Design Functions

```matlab
% NTF Synthesis
[ntf, stf] = synthesizeNTF(order, OSR, opt, H_inf, f0);

% NTF Realization
[a, g, b, c] = realizeNTF(ntf, form);

% ABCD Construction
ABCD = stuffABCD(a, g, b, c, form);

% Simulate (toolbox function)
[v, xn, xmax, y] = simulateDSM(u, ABCD, nlev, x0);
```

### Analysis Functions

```matlab
% Plot NTF
plotNTF(ntf, OSR, f0);

% Predict SNR
snr = predictSNR(ntf, OSR);

% Calculate noise gain
rms = rmsGain(ntf, f0, fB);

% Find stable input range
[amp, snr] = simulateSNR(ntf, OSR, [], [], nlev);
```

### Plotting

Reference designs create visible figure windows. To show plots:
```matlab
% Create visible figure (default behavior)
figure('Name', 'My Plot');  % or simply: figure;

% For headless/batch mode, use:
figure('Visible', 'off');
```

When running interactively in MATLAB, figures will display automatically.

## Example Specifications

### Example 1: 2nd-Order CIFF

```matlab
order = 2;
OSR = 64;
H_inf = 1.5;
opt = 1;
fs = 1e6;
n_bits = 5;
A_in = 0.5;

% Expected: SNR ~70 dB
```

### Example 2: 3rd-Order CIFF

```matlab
order = 3;
OSR = 16;
H_inf = 2.0;
opt = 1;
fs = 2e6;
n_bits = 5;
A_in = 0.2;

% Expected: SNR ~74 dB
```

### Example 3: 4th-Order CIFF (Aggressive)

```matlab
order = 4;
OSR = 16;
H_inf = 4.0;  % Note: aggressive, needs verification
opt = 1;
fs = 2e6;
n_bits = 5;
A_in = 0.1;   % Lower input for stability

% Expected: SNR ~90 dB
```

## SNR Calculation (Critical!)

**CORRECT:** Find signal in **OUTPUT spectrum**, **limit search to in-band only**
```matlab
V_out = fft(v .* w) / (N/4);
fB_bins = ceil(N / (2*OSR));
[~, sig_idx] = max(V_out_mag(2:fB_bins));  % Search only in-band
sig_power = sum(V_out_mag(sig_bins).^2);
noise_power = sum(V_out_mag(noise_bins).^2);
SNR = 10*log10(sig_power / noise_power);
```

**WRONG:** Searching beyond signal band
```matlab
[~, sig_idx] = max(V_out_mag(2:N/2));  % May find spur outside band!
```

**WRONG:** Finding signal in noise spectrum
```matlab
noise = v - u;
V_noise = fft(noise .* w);  % DON'T find signal here!
```

### DC Leakage from Hann Window

When using Hann windowing, DC offset (from DAC mismatch) leaks into bin 2.
**Exclude bin 2 from noise calculation:**

```matlab
% Noise calculation excludes DC (bin 1) and low-freq leakage (bin 2)
fB_bins = ceil(N / (2*OSR));
noise_bins = setdiff(3:fB_bins, [sig_bins, harmonic_bins]);  % Start from bin 3
```

### SNR vs SNDR

- **SNR:** Signal / (Noise only). Excludes harmonic distortion.
  - Signal search: limited to in-band bins (2:fB_bins)
  - Noise bins: Everything in-band except DC, signal, harmonics, and bin 2
- **SNDR:** Signal / (Noise + Distortion). Includes harmonics.
  - Signal search: limited to in-band bins (2:fB_bins)
  - Noise bins: Everything in-band except DC, signal, and bin 2

## Troubleshooting

### "synthesizeNTF not found"
Add Delta Sigma Toolbox to path (relative to your script):
```matlab
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));
```

Or for global installation:
```matlab
addpath('/usr/local/MATLAB/toolbox/deltasig');
```

### Simulation unstable (states diverge)
- Reduce input amplitude (try 0.1V instead of 0.5V)
- Check H_inf not too high for order (Lee's rule: < 2 for stability)
- Verify coefficients are reasonable (a values ~1, not 1e10)

### Poor SNR
- Check signal is in passband (not at harmonic)
- Verify OSR calculation: fB = fs/(2*OSR)
- Ensure windowing before FFT

### Plots not showing
Reference designs create visible figures by default. If plots don't appear:
- **Linux**: Ensure `$DISPLAY` is set or use `matlab -desktop`
- **SSH**: Enable X11 forwarding: `ssh -X user@host`
- **Headless**: Plots are also saved as PNG files in the working directory

## Reference Designs

Located in `references/` directory. All scripts include automatic toolbox path setup:

```matlab
% At the top of each script
addpath(fullfile(fileparts(mfilename('fullpath')), 'dstoolbox'));
```

**Folder Structure:**
```
references/
├── dstoolbox/                 # Delta Sigma Toolbox
├── design_4th_order_ciff.m    # Full 4th-order example (plots visible)
├── dsm_4th_order_simple.m     # Headless version (no plots)
└── dsm_quick_design.m         # Quick-start template
```

### 4th-Order CIFF Design (2MSPS, OSR=16, 5-bit, H_inf=4.0)

See `references/design_4th_order_ciff.m` for complete example.

**Specifications:**
| Parameter | Value |
|-----------|-------|
| Order | 4 |
| fs | 2 MHz |
| OSR | 16 (BW = 62.5 kHz) |
| H_inf | 4.0 |
| Topology | CIFF |
| Quantizer | 5-bit (32 levels) |

**Expected Results:**
| Input | SNR | ENOB |
|-------|-----|------|
| 0.1V | ~90 dB | ~14.6 bits |
| 0.5V | ~104 dB | ~17.0 bits |

## ADC and DAC Component Modeling

In addition to complete DSM designs, this skill supports component-level modeling of flash ADC quantizers and thermometer DACs.

### Flash ADC Quantizer with Thermometer Code Output

A flash ADC uses a bank of comparators operating in parallel. The output is thermometer code (contiguous 1s followed by 0s), which can be converted to binary output.

```matlab
%% Flash ADC Quantizer
% Parameters
v_fs = 1.0;   % +/- 1V full scale
n_bits = 3;   % 3-bit resolution

% Input signal
v_in = linspace(-1, 1, 100);

% Quantize
[thermometer, binary, thresholds] = flash_adc_quantizer(v_in, v_fs, n_bits);

% thermometer: [samples x (2^n_bits - 1)] matrix of 0/1
% binary: quantized output voltages (row vector)
% thresholds: comparator threshold voltages
```

**Thermometer Code:** For a 3-bit flash ADC (7 comparators):
- `0000000` → Level 0 → -1.00V
- `0001111` → Level 4 → +0.14V  
- `1111111` → Level 7 → +1.00V

Each bit represents a comparator: `1` if `vin > threshold`, `0` otherwise.

### Thermometer DAC

Converts thermometer code back to analog voltage. Each bit has weight = `2*v_fs / (2^n_bits - 1)`.

```matlab
%% Thermometer DAC
% Convert thermometer code back to analog
v_dac = thermometer_dac(thermometer, v_fs, n_bits);

% v_dac matches the binary output from the ADC
% Both are row vectors [1 x samples]
```

**DAC Operation:**
- Bit = 1: contributes `+weight`
- Bit = 0: contributes `-weight`
- Output = (sum of contributions) / 2 → scaled to ±v_fs

### ADC-DAC Chain Verification

```matlab
%% Complete ADC-DAC chain
% Flash ADC
[therm, binary_adc, ~] = flash_adc_quantizer(v_in, 1.0, 3);

% Thermometer DAC
binary_dac = thermometer_dac(therm, 1.0, 3);

% Verify they match
max_error = max(abs(binary_adc - binary_dac));  % Should be ~0
```

### DAC Unit Cell Mismatch Modeling

Real-world DACs suffer from unit cell mismatch due to manufacturing variations. This causes non-uniform step sizes (DNL errors) and deviation from the ideal transfer curve (INL errors).

**Mismatch Model:**
```
Actual weight = nominal_weight + error
where error ~ N(0, (mismatch_pct * nominal_weight)^2)
```

```matlab
%% DAC with Unit Cell Mismatch
% Parameters
mismatch_pct = 0.02;  % 2% standard deviation
seed = 42;            % Random seed for reproducibility

% Run ADC
[thermometer, binary_ideal, ~] = flash_adc_quantizer(v_in, 1.0, 3);

% Run DAC with mismatch
[v_mismatch, bit_weights] = thermometer_dac_mismatch(thermometer, 1.0, 3, mismatch_pct, seed);

% bit_weights contains actual weights for each unit cell
```

**Mismatch Effects:**
- **DNL (Differential Nonlinearity)**: Deviation of step size from 1 LSB
- **INL (Integral Nonlinearity)**: Deviation of transfer curve from ideal
- **Harmonic Distortion**: Creates spurs in output spectrum
- **Reduced SFDR**: Limits dynamic range of multi-bit DACs

**Example Output with 2% Mismatch:**
```
Nominal weight: 0.2857V
Cell 1: 0.2912V (+1.9% deviation)
Cell 2: 0.2789V (-2.4% deviation)
...
INL (max): 0.15 LSB
DNL (max): 0.08 LSB
```

**Available Functions:**
- `thermometer_dac_mismatch.m` - DAC with unit cell mismatch
- `thermometer_dac_mismatch_demo.m` - Demonstration of mismatch effects with INL/DNL analysis

### Dynamic Weighted Averaging (DWA)

DWA is a mismatch-shaping technique that rotates which unit cells are selected for a given thermometer code, averaging out mismatch errors over time.

**Algorithm:**
- For each sample with `k` ones in the thermometer code:
  1. Select `k` consecutive cells starting from the current position
  2. Wrap around if needed (circular buffer)
  3. Next start position = position after last selected cell

**Example with 7 cells (3-bit), code = 3:**
```
Sample 1: cells 1,2,3 activated → next start = 4
Sample 2: cells 4,5,6 activated → next start = 7
Sample 3: cells 7,1,2 activated (wrap) → next start = 3
Sample 4: cells 3,4,5 activated → next start = 6
```

```matlab
%% DWA DAC with Mismatch
[v_out, start_idx] = thermometer_dac_dwa(thermometer, v_fs, n_bits, bit_weights, start_idx);
```

**Benefits:**
- Equalizes cell usage across all unit cells
- Shapes mismatch noise out of signal band
- Pushes distortion to higher frequencies (attenuated by NTF)
- Typical improvement: 10-20 dB in SNDR with 0.5-2% mismatch

**Implementation in DSM:**
See `design_4th_order_ciff_with_dac_mismatch.m` for complete example comparing:
- Ideal DAC (no mismatch)
- Static DAC with mismatch
- DWA DAC with same mismatch

**Key Implementation Detail:**
The DWA DAC creates a **rotated** thermometer pattern for feedback, while the ADC output uses the original pattern:
```matlab
% ADC output uses original thermometer
them = flash_adc_quantizer(y, V_fs, n_bits);
v_adc(i) = binary_output;

% DAC feedback uses rotated pattern
dwa_therm = rotate_thermometer(therm, dwa_start_idx);
v_dac(i) = calculate_dac_output(dwa_therm, bit_weights);
```

**Available Functions:**
- `thermometer_dac_dwa.m` - DWA DAC implementation
- `thermometer_dac_dwa_demo.m` - Standalone DWA demo with SNDR comparison

## References

1. Pavan, Schreier, Temes - "Understanding Delta-Sigma Data Converters" (2nd Ed.), Appendix B
2. Schreier, R. - "The Delta-Sigma Toolbox" (MATLAB Central)
3. Delta Sigma Toolbox: https://www.mathworks.com/matlabcentral/fileexchange/19-delta-sigma-toolbox

---
*Native MATLAB workflow - no Python dependencies*
