"""
Plot coefficient ROC diagnostic from diag_coeff_evo.csv.

Two figures are produced:

  diag_coeff_evo_ratio.png  – |C_{m+1}/C_m| vs m, one line per ψ.
                               Horizontal reference lines: e (ae bound) and
                               ae²/ρ_evo(ψ) (evolute bound) for each angle.

  diag_coeff_evo_abs.png    – log|C_m| vs m, one line per ψ.
                               Slope = log(1/R) confirms the ROC.

Usage:
    python python/plot_diag_coeff_evo.py
"""

import csv
import io
import math
import os
import re
import numpy as np
import matplotlib.pyplot as plt

CSV_PATH = os.path.join(os.path.dirname(__file__), "..", "test_data", "diag_coeff_evo.csv")
OUT_DIR  = os.path.join(os.path.dirname(__file__), "..", "test_data")

ABS_FLOOR = 1e-300

# ---------------------------------------------------------------------------
# Load CSV
# ---------------------------------------------------------------------------
with open(CSV_PATH, newline="") as f:
    raw = f.read()

meta_match = re.search(r"#\s*a=([\d.eE+\-]+)\s+b=([\d.eE+\-]+)", raw)
if not meta_match:
    raise ValueError(f"No '# a=... b=...' metadata line found in {CSV_PATH}")
ELLIPSE_A = float(meta_match.group(1))
ELLIPSE_B = float(meta_match.group(2))

e2   = 1.0 - (ELLIPSE_B / ELLIPSE_A) ** 2
e    = math.sqrt(e2)
ae   = ELLIPSE_A * e
ae2  = ELLIPSE_A * e2

csv_lines = (l for l in raw.splitlines(keepends=True) if not l.startswith("#"))
reader = csv.DictReader(io.StringIO("".join(csv_lines)))

# data[psi_deg] -> list of (m, C_m, abs_C_m, ratio) sorted by m
data = {}
for row in reader:
    psi   = float(row["psi_deg"])
    m     = int(row["m"])
    cm    = float(row["C_m"])
    acm   = float(row["abs_C_m"])
    ratio = float(row["ratio"])
    data.setdefault(psi, []).append((m, cm, acm, ratio))

for psi in data:
    data[psi].sort(key=lambda t: t[0])

psi_list = sorted(data.keys())


def evolute_rho(psi_deg):
    """Evolute radius at psi_deg (same units as a)."""
    psi  = math.radians(psi_deg)
    ac   = abs(ELLIPSE_A * math.cos(psi))
    bs   = abs(ELLIPSE_B * math.sin(psi))
    s    = ac ** (2/3) + bs ** (2/3)
    return (ELLIPSE_A**2 - ELLIPSE_B**2) / s ** 1.5


# ---------------------------------------------------------------------------
# Figure 1: ratio |C_{m+1}/C_m| vs m
# ---------------------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(10, 6))
fig1.suptitle(
    f"Coefficient ratio |C_{{m+1}}/C_m| vs m\n"
    f"a={ELLIPSE_A}, b/a={ELLIPSE_B/ELLIPSE_A:.4g},  "
    f"e={e:.4f},  ae={ae:.4f},  ae²={ae2:.4f}"
)

colors = plt.cm.tab10(np.linspace(0, 1, len(psi_list)))

for psi, color in zip(psi_list, colors):
    rows  = data[psi]
    ms    = [r[0] for r in rows if r[3] >= 0]
    rats  = [r[3] for r in rows if r[3] >= 0]
    ax1.plot(ms, rats, color=color, label=f"ψ={psi:.4g}°")

    # Reference: ae bound (ratio → e) and evolute bound (ratio → ae²/ρ_evo)
    rho_evo = evolute_rho(psi)
    evo_ref = ae2 / rho_evo
    psi_rad = math.radians(psi)
    # Only draw ae reference if ae < rho_evo (ae-dominated)
    binding = "ae" if ae < rho_evo else "evo"
    if binding == "ae":
        ax1.axhline(e, color=color, linestyle=":", linewidth=0.8)
    else:
        ax1.axhline(evo_ref, color=color, linestyle=":", linewidth=0.8)

# Permanent label lines
ax1.axhline(e, color="gray", linestyle="--", linewidth=1.0, label=f"e={e:.4f}  (ae bound)")

ax1.set_xlabel("m")
ax1.set_ylabel("|C_m| / |C_{m-1}|")
ax1.set_ylim(0, 1.5)
ax1.legend(fontsize=8)
ax1.grid(True, alpha=0.3)

out1 = os.path.join(OUT_DIR, "diag_coeff_evo_ratio.png")
fig1.savefig(out1, dpi=150, bbox_inches="tight")
print(f"Saved {out1}")
plt.close(fig1)

# ---------------------------------------------------------------------------
# Figure 2: log|C_m| vs m
# ---------------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(10, 6))
fig2.suptitle(
    f"log|C_m| vs m  (slope = log(1/R))\n"
    f"a={ELLIPSE_A}, b/a={ELLIPSE_B/ELLIPSE_A:.4g}"
)

for psi, color in zip(psi_list, colors):
    rows  = data[psi]
    ms    = [r[0] for r in rows]
    acms  = [max(r[2], ABS_FLOOR) for r in rows]
    ax2.semilogy(ms, acms, color=color, label=f"ψ={psi:.4g}°")

ax2.set_xlabel("m")
ax2.set_ylabel("|C_m|")
ax2.legend(fontsize=8)
ax2.grid(True, which="both", alpha=0.3)

out2 = os.path.join(OUT_DIR, "diag_coeff_evo_abs.png")
fig2.savefig(out2, dpi=150, bbox_inches="tight")
print(f"Saved {out2}")
plt.close(fig2)
