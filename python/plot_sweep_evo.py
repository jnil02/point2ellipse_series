"""
Plot convergence results from sweep_evo.csv.

Two subplots for a chosen psi angle:
  Left:  |phi error| vs rho  for several truncation orders N
  Right: |phi error| vs N    for several rho values (inside and outside evolute)

The evolute boundary rho_evo is marked on both plots.

Usage:
    python python/plot_sweep_evo.py
"""

import csv
import math
import os
import matplotlib.pyplot as plt

CSV_PATH = os.path.join(os.path.dirname(__file__), "..", "test_data", "sweep_evo.csv")
OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "test_data", "sweep_evo.png")

PSI_PLOT  = 45  # which psi angle (degrees) to inspect
N_COUNT   = 6   # how many N curves in the left plot  (evenly spaced from data)
RHO_COUNT = 5   # how many rho curves in the right plot (evenly spaced from data)

# ---------------------------------------------------------------------------
# Load CSV into a nested dict: data[psi_deg][N] = list of (rho, rho_evo, phi_err, h_err)
# ---------------------------------------------------------------------------
data = {}   # data[psi_deg][N] -> list of (rho, rho_evo, phi_err, h_err)

with open(CSV_PATH, newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        psi  = int(row["psi_deg"])
        N    = int(row["N"])
        rho      = float(row["rho"])
        rho_evo  = float(row["rho_evo"])
        phi_err  = float(row["phi_err"])
        h_err    = float(row["h_err"])
        data.setdefault(psi, {}).setdefault(N, []).append((rho, rho_evo, phi_err, h_err))

# Sort each list by rho.
for psi in data:
    for N in data[psi]:
        data[psi][N].sort(key=lambda t: t[0])

# ---------------------------------------------------------------------------
# Extract data for the chosen psi.
# ---------------------------------------------------------------------------
if PSI_PLOT not in data:
    raise ValueError(f"psi={PSI_PLOT}° not found in {CSV_PATH}. "
                     f"Available: {sorted(data.keys())}")

psi_data = data[PSI_PLOT]
N_max    = max(psi_data.keys())

# rho_evo is the same for all N at a given psi — take it from N=1.
rho_evo = psi_data[1][0][1]

# Unique rho values (from any N).
rhos_all = [t[0] for t in psi_data[1]]

# Pick N_COUNT evenly spaced N values from whatever the data contains.
Ns_all = sorted(psi_data.keys())
def _pick(lst, count):
    n = len(lst)
    if n <= count:
        return lst[:]
    return [lst[round(i * (n - 1) / (count - 1))] for i in range(count)]

N_LINES  = _pick(Ns_all, N_COUNT)
rho_pick = _pick(rhos_all, RHO_COUNT)

# ---------------------------------------------------------------------------
# Plot.
# ---------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
fig.suptitle(
    f"Inside-evolute series convergence  "
    f"(a=1, b/a=0.5,  ψ={PSI_PLOT}°,  ρ_evo={rho_evo:.4f})"
)

# --- Left: phi_err vs rho for selected N ---
for N in N_LINES:
    if N not in psi_data:
        continue
    rhos     = [t[0] for t in psi_data[N]]
    phi_errs = [max(t[2], 1e-50) for t in psi_data[N]]   # guard log(0)
    ax1.semilogy(rhos, phi_errs, label=f"N={N}")

ax1.axvline(rho_evo, color="k", linestyle="--", linewidth=1, label=f"ρ_evo={rho_evo:.3f}")
ax1.set_xlabel("ρ")
ax1.set_ylabel("|φ error| [rad]")
ax1.set_title("Error vs ρ")
ax1.legend()
ax1.grid(True, which="both", alpha=0.3)

for rho_val in rho_pick:
    phi_errs = []
    for N in Ns_all:
        # Find the row closest to rho_val.
        row = min(psi_data[N], key=lambda t: abs(t[0] - rho_val))
        phi_errs.append(max(row[2], 1e-50))
    inside = "in" if rho_val < rho_evo else "out"
    ax2.semilogy(Ns_all, phi_errs, marker="o", label=f"ρ={rho_val:.3f} ({inside})")

ax2.set_xlabel("N  (= K)")
ax2.set_ylabel("|φ error| [rad]")
ax2.set_title("Error vs N")
ax2.legend()
ax2.grid(True, which="both", alpha=0.3)

plt.tight_layout()
plt.savefig(OUT_PATH, dpi=150)
print(f"Saved {OUT_PATH}")
plt.show()
