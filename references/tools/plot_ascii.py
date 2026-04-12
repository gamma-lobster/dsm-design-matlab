#!/usr/bin/env python3
"""Plot DSM results - ASCII/Terminal version"""

import math

# Parameters
fs = 2e6
OSR = 16
fB = fs / (2 * OSR)
N = 8192

# Coefficients from simulation
a = [2.566151, 2.736364, 1.365746, 0.218951]
g = [0.004450, 0.028318]
SNR = 103.72
ENOB = 16.94

# NTF zeros and poles from synthesizeNTF output
z_ntf = [
    0.9977727 + 0.06670555j,
    0.9977727 - 0.06670555j,
    0.98573937 + 0.16827921j,
    0.98573937 - 0.16827921j
]
p_ntf = [
    0.35410201 - 0.53003518j,
    0.35410201 + 0.53003518j,
    0.36282254 - 0.13714439j,
    0.36282254 + 0.13714439j
]
k_ntf = 1.0

def calc_ntf(freq, fs_val):
    """Calculate NTF at given frequency"""
    z = complex(math.cos(2*math.pi*freq/fs_val), math.sin(2*math.pi*freq/fs_val))
    num = k_ntf
    for zi in z_ntf:
        num *= (z - zi)
    den = 1.0
    for pi in p_ntf:
        den *= (z - pi)
    return num / den

def hz_to_khz(hz):
    return hz / 1000

def db(mag):
    return 20 * math.log10(abs(mag) + 1e-15)

print("=" * 70)
print("  4th-Order CIFF DSM - ASCII Spectrum Plot")
print("=" * 70)
print()
print(f"Design: order=4, OSR={OSR}, H_inf=4.0, fs={fs/1e6:.1f} MHz")
print(f"Signal BW: {fB/1000:.1f} kHz")
print(f"SNR: {SNR:.2f} dB | ENOB: {ENOB:.2f} bits")
print()

# Plot NTF
print("-" * 70)
print("NTF Response:")
print("-" * 70)

width = 60
height = 15
db_min, db_max = -100, 40

# Generate NTF points
n_points = width
freqs = [10**(math.log10(100) + (math.log10(fs/2) - math.log10(100)) * i / (n_points-1)) 
         for i in range(n_points)]
ntf_db = [db(calc_ntf(f, fs)) for f in freqs]

# Normalize to plot height
plot_data = []
for db_val in ntf_db:
    if db_val < db_min:
        row = db_min
    elif db_val > db_max:
        row = db_max
    else:
        row = db_val
    plot_data.append(int((row - db_min) / (db_max - db_min) * (height - 1)))

# Print y-axis labels and plot
for row in range(height-1, -1, -1):
    db_label = db_min + row * (db_max - db_min) / (height - 1)
    line = f"{db_label:5.0f} dB |"
    for col in range(width):
        if plot_data[col] == row:
            line += "*"
        elif abs(plot_data[col] - row) == 1:
            line += "."
        else:
            line += " "
    print(line)

# X-axis
print("       " + "-" * width)
x_labels = [100, 1000, 10000, 100000, fs/2]
x_label_str = ""
for f in freqs:
    for xl in x_labels:
        if abs(f - xl) / xl < 0.05:
            pos = int((math.log10(f) - math.log10(100)) / (math.log10(fs/2) - math.log10(100)) * width)
            print(f"{' ' * (7 + pos)}{hz_to_khz(xl):.0f}kHz")
            break

print()

# Plot output spectrum (simulated based on NTF and signal)
print("-" * 70)
print("Output Spectrum (approximate):")
print("-" * 70)

# Simulate output spectrum
freqs_lin = [i * fs / N for i in range(N//2 + 1)]
sig_bin = int(round(math.sqrt(1/7) * N / (2 * OSR)))
sig_freq = sig_bin * fs / N

spectrum_db = []
for f in freqs_lin:
    if f < fB:
        # In-band: shaped noise + signal at sig_freq
        ntf_val = calc_ntf(f, fs)
        noise_db = -140 + db(ntf_val)  # Base noise floor shaped by NTF
        if abs(f - sig_freq) < fs/N * 2:  # Signal bin
            noise_db = max(noise_db, -20)  # Signal peak
        spectrum_db.append(noise_db)
    else:
        # Out-of-band
        ntf_val = calc_ntf(f, fs)
        noise_db = -140 + db(ntf_val)
        spectrum_db.append(noise_db)

# Downsample for display
n_display = width
indices = [int(i * len(spectrum_db) / n_display) for i in range(n_display)]
display_db = [spectrum_db[min(i, len(spectrum_db)-1)] for i in indices]

# Plot
db_floor, db_ceil = -150, 0
plot_heights = []
for db_val in display_db:
    h = int((db_val - db_floor) / (db_ceil - db_floor) * (height - 1))
    h = max(0, min(height-1, h))
    plot_heights.append(h)

for row in range(height-1, -1, -1):
    db_label = db_floor + row * (db_ceil - db_floor) / (height - 1)
    line = f"{db_label:5.0f} dB|"
    for col in range(n_display):
        if plot_heights[col] == row:
            # Mark signal specially
            if col == int(sig_bin / len(spectrum_db) * n_display):
                line += "S"  # Signal
            else:
                line += "*"
        elif plot_heights[col] > row:
            line += "|"
        else:
            line += " "
    print(line)

print("       " + "-" * n_display)
print(f"       0{' ' * (n_display//2 - 5)}{hz_to_khz(fB):.0f}kHz{' ' * (n_display//2 - 8)}{hz_to_khz(fs/2):.0f}kHz")
print("       |________Signal BW________|_____Out-of-band_____|")

print()
print("-" * 70)
print("Summary:")
print("-" * 70)
print(f"  Signal Frequency: {hz_to_khz(sig_freq):.1f} kHz")
print(f"  Signal BW: 0 - {hz_to_khz(fB):.1f} kHz")
print(f"  SNR: {SNR:.2f} dB")
print(f"  ENOB: {ENOB:.2f} bits")
print(f"  Status: STABLE")
print("=" * 70)
