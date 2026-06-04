"""
Plot convergence results from sweep_evo_m.csv.

For each error column in the CSV (any column that is not psi_deg / rho /
rho_evo / N) two figures are produced:

  <col>.png        – two subplots for a chosen psi angle:
                       Left:  |error| vs rho  for several truncation orders N
                       Right: |error| vs N    for several rho values
  <col>_polar.png  – polar heat map of the decay rate (d log10|err| / dN)
                     across the full (psi, rho) grid with the evolute overlaid.

Usage:
    python python/plot_sweep_evo.py
"""

import csv
import io
import math
import os
import re
import numpy as np
import matplotlib.pyplot as plt

CSV_PATH = os.path.join(os.path.dirname(__file__), "..", "test_data", "sweep_evo_m.csv")
OUT_DIR  = os.path.join(os.path.dirname(__file__), "..", "test_data")

PSI_PLOT  = 80  # which psi angle (degrees) to inspect
N_COUNT   = 6   # how many N curves in the left plot  (evenly spaced from data)
RHO_COUNT = 5   # how many rho curves in the right plot (evenly spaced from data)
TAIL      = 5   # high-N points used for slope estimate in the heat map
ERR_FLOOR = 1e-300  # log-plot guard; set to match ~10^-(BITS*log10(2)) for chosen precision

NON_ERR_COLS = {"psi_deg", "rho", "rho_evo", "N"}

# ---------------------------------------------------------------------------
# Load CSV — first line is "# a=... b=..." metadata written by sweep_evo.cpp
# ---------------------------------------------------------------------------
with open(CSV_PATH, newline="") as f:
    raw = f.read()

meta_match = re.search(r"#\s*a=([\d.eE+\-]+)\s+b=([\d.eE+\-]+)", raw)
if not meta_match:
    raise ValueError(f"No '# a=... b=...' metadata line found in {CSV_PATH}")
ELLIPSE_A = float(meta_match.group(1))
ELLIPSE_B = float(meta_match.group(2))

csv_lines = (l for l in raw.splitlines(keepends=True) if not l.startswith("#"))
reader = csv.DictReader(io.StringIO("".join(csv_lines)))

err_cols = [c for c in reader.fieldnames if c not in NON_ERR_COLS]
if not err_cols:
    raise ValueError("No error columns found in CSV.")

# data[psi_deg][N] -> list of (rho, rho_evo, {col: value, ...})
data = {}
for row in reader:
    psi     = float(row["psi_deg"])
    N       = int(row["N"])
    rho     = float(row["rho"])
    rho_evo = float(row["rho_evo"])
    errs    = {col: float(row[col]) for col in err_cols}
    data.setdefault(psi, {}).setdefault(N, []).append((rho, rho_evo, errs))

for psi in data:
    for N in data[psi]:
        data[psi][N].sort(key=lambda t: t[0])

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _pick(lst, count):
    """Pick `count` evenly-spaced elements from lst."""
    n = len(lst)
    if n <= count:
        return lst[:]
    return [lst[round(i * (n - 1) / (count - 1))] for i in range(count)]

def _edges(a):
    """Cell edges for pcolormesh: clamp left, extend right by half step."""
    d = np.diff(a)
    mid = (a[:-1] + a[1:]) / 2
    return np.r_[a[0], mid, a[-1] + d[-1] / 2]

# ---------------------------------------------------------------------------
# Per-column plots
# ---------------------------------------------------------------------------
psi_degs_all = sorted(data.keys())
Ns_heat      = sorted(data[psi_degs_all[0]].keys())
rhos_heat    = [t[0] for t in data[psi_degs_all[0]][Ns_heat[0]]]
ae_val       = ELLIPSE_A * math.sqrt(1 - (ELLIPSE_B / ELLIPSE_A) ** 2)

psi_plot_actual = min(data.keys(), key=lambda p: abs(p - PSI_PLOT))
if abs(psi_plot_actual - PSI_PLOT) > 5:
    raise ValueError(f"No psi near {PSI_PLOT}° in {CSV_PATH}. "
                     f"Available: {sorted(data.keys())}")

psi_data = data[psi_plot_actual]
Ns_all   = sorted(psi_data.keys())
rho_evo  = psi_data[1][0][1]
rhos_all = [t[0] for t in psi_data[1]]
N_LINES  = _pick(Ns_all, N_COUNT)
rho_pick = _pick(rhos_all, RHO_COUNT)

for col in err_cols:
    # --- Cartesian figure ---
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        f"Inside-evolute series convergence — {col}\n"
        f"a={ELLIPSE_A}, b/a={ELLIPSE_B/ELLIPSE_A:.4g},  "
        f"ψ={psi_plot_actual:.4g}°,  ρ_evo={rho_evo:.4f}"
    )

    for N in N_LINES:
        if N not in psi_data:
            continue
        rhos = [t[0] for t in psi_data[N]]
        errs = [max(t[2][col], ERR_FLOOR) for t in psi_data[N]]
        ax1.semilogy(rhos, errs, label=f"N={N}")

    ax1.axvline(rho_evo, color="k", linestyle="--", linewidth=1,
                label=f"ρ_evo={rho_evo:.3f}")
    ax1.set_xlabel("ρ")
    ax1.set_ylabel(f"|{col}|")
    ax1.set_title("Error vs ρ")
    ax1.legend()
    ax1.grid(True, which="both", alpha=0.3)

    for rho_val in rho_pick:
        errs = []
        for N in Ns_all:
            row = min(psi_data[N], key=lambda t: abs(t[0] - rho_val))
            errs.append(max(row[2][col], ERR_FLOOR))
        inside = "in" if rho_val < rho_evo else "out"
        ax2.semilogy(Ns_all, errs, marker="o", label=f"ρ={rho_val:.3f} ({inside})")

    ax2.set_xlabel("N  (= K)")
    ax2.set_ylabel(f"|{col}|")
    ax2.set_title("Error vs N")
    ax2.legend()
    ax2.grid(True, which="both", alpha=0.3)

    plt.tight_layout()
    out_path = os.path.join(OUT_DIR, f"sweep_evo_{col}.png")
    plt.savefig(out_path, dpi=150)
    print(f"Saved {out_path}")
    plt.close(fig)

    # --- Polar heat-map figure ---
    rate = np.full((len(psi_degs_all), len(rhos_heat)), np.nan)

    for i, psi in enumerate(psi_degs_all):
        Ns_psi = sorted(data[psi].keys())
        for j, rho_val in enumerate(rhos_heat):
            log_errs, ns_used = [], []
            for N in Ns_psi:
                row = min(data[psi][N], key=lambda t: abs(t[0] - rho_val))
                e = row[2][col]
                if e > ERR_FLOOR:
                    log_errs.append(math.log10(e))
                    ns_used.append(float(N))
            if len(log_errs) >= 2:
                rate[i, j] = np.polyfit(ns_used[-TAIL:], log_errs[-TAIL:], 1)[0]

    thetas   = np.radians(np.array(psi_degs_all))
    rhos_arr = np.array(rhos_heat)

    fig_polar, ax_polar = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 6))
    fig_polar.suptitle(
        f"Decay rate — {col}  (d log₁₀|err| / dN)\n"
        "blue = fast convergence,  red = divergence"
    )

    vmax = np.nanmax(np.abs(rate))
    mesh = ax_polar.pcolormesh(
        _edges(thetas), _edges(rhos_arr), rate.T,
        cmap="RdBu_r", vmin=-vmax, vmax=vmax, shading="flat"
    )
    fig_polar.colorbar(mesh, ax=ax_polar, pad=0.1,
                       label=f"d log₁₀|{col}| / dN")

    T, R = np.meshgrid(thetas, rhos_arr)
    ax_polar.contour(T, R, rate.T, levels=15, colors="k", linewidths=0.4, alpha=0.5)
    cs0 = ax_polar.contour(T, R, rate.T, levels=[0], colors="k", linewidths=0.4)
    ax_polar.clabel(cs0, fmt="0", fontsize=8)

    evo_rhos = [data[psi][Ns_heat[0]][0][1] for psi in psi_degs_all]
    ax_polar.plot(thetas, evo_rhos, "k--", linewidth=1.5, label="evolute")

    psi_all_rad = np.linspace(0, np.pi / 2, 200)
    ax_polar.plot(psi_all_rad, np.full_like(psi_all_rad, ae_val),
                  "b--", linewidth=1.5, label=f"ae = {ae_val:.3f}")
    ax_polar.legend(loc="upper left", bbox_to_anchor=(1.15, 1.05))

    out_path_polar = os.path.join(OUT_DIR, f"sweep_evo_{col}_polar.png")
    fig_polar.savefig(out_path_polar, dpi=150, bbox_inches="tight")
    print(f"Saved {out_path_polar}")
    plt.close(fig_polar)
