"""Mixing-layer profiles and the paper's key numbers (Hillier & Arregui 2019, Fig. 6).

Usage (from the repo root, with .venv active):
    python scripts/profiles.py paper2d/output/kh_hlld 30
    python scripts/profiles.py paper2d/output/kh_hlld 20 30   # time-average snapshots 20..30

Computes x-averaged profiles <rho>(y), <v_x>(y) and the mean-flow kinetic energy
1/2 <rho><v_x>^2, finds the mixing-layer edges, and compares with the model:
  - mean density in the layer        -> sqrt(rho1 rho2)
  - layer asymmetry |y_low| / y_up   -> sqrt(rho2/rho1)   (shift to the light side)
  - fraction of mean-flow KE removed -> 0.53 (model), 0.58 (paper's 3D simulation)
Saves profiles.png.
"""
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_temperature import load_uniform  # noqa: E402

RHO1, RHO2 = 1.0, 10.0   # light (y<0) and dense (y>0) densities
EDGE_TOL = 0.01          # paper: edge where <rho> departs from background by 1%


def mean_profiles(fname):
    """x-averaged density and momentum-weighted x-velocity."""
    x, y, f, t = load_uniform(fname)
    rho = f["rho"].mean(axis=0)
    mom = f["m1"].mean(axis=0)
    return y, rho, mom, t


def layer_edges(y, rho):
    """Lowest y where <rho> rises 1% above rho1, highest y where it is 1% below rho2."""
    inside = (rho > RHO1 * (1 + EDGE_TOL)) & (rho < RHO2 * (1 - EDGE_TOL))
    if not inside.any():
        raise RuntimeError("no mixing layer found")
    return y[inside].min(), y[inside].max()


def band_kinetic_energy(y, rho, mom, lo, hi):
    """Mean-flow KE in [lo, hi], measured in that band's zero-momentum frame."""
    band = (y >= lo) & (y <= hi)
    dy = y[1] - y[0]
    v = mom[band] / rho[band]
    v_frame = mom[band].sum() / rho[band].sum()   # band moves with this velocity
    return 0.5 * np.sum(rho[band] * (v - v_frame) ** 2) * dy, v_frame


def main():
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    base = sys.argv[1]
    n_first = int(sys.argv[2])
    n_last = int(sys.argv[3]) if len(sys.argv) == 4 else n_first

    # initial state, for the KE reference
    y0, rho0, mom0, _ = mean_profiles(f"{base}0000.dat")

    # (time-)averaged profiles
    rho_sum = mom_sum = None
    times = []
    for n in range(n_first, n_last + 1):
        y, rho, mom, t = mean_profiles(f"{base}{n:04d}.dat")
        rho_sum = rho if rho_sum is None else rho_sum + rho
        mom_sum = mom if mom_sum is None else mom_sum + mom
        times.append(t)
    nsnap = n_last - n_first + 1
    rho, mom = rho_sum / nsnap, mom_sum / nsnap
    vx = mom / rho

    lo, hi = layer_edges(y, rho)
    band = (y >= lo) & (y <= hi)
    rho_mean_layer = rho[band].mean()

    # same band at t = 0 (interpolate in case the grid level differs)
    rho0_i = np.interp(y, y0, rho0)
    mom0_i = np.interp(y, y0, mom0)
    ke0, _ = band_kinetic_energy(y, rho0_i, mom0_i, lo, hi)
    ke, v_frame = band_kinetic_energy(y, rho, mom, lo, hi)

    # model predictions
    rho_model = np.sqrt(RHO1 * RHO2)
    asym_model = np.sqrt(RHO2 / RHO1)

    tlabel = f"t = {times[0]:.0f}" if nsnap == 1 else f"t = {times[0]:.0f}-{times[-1]:.0f}"
    print(f"Snapshot(s): {tlabel}")
    print(f"Layer edges:            y_low = {lo:+.3f},  y_up = {hi:+.3f}  (width {hi - lo:.3f})")
    print(f"Asymmetry |y_low|/y_up: {abs(lo) / hi:6.2f}   model sqrt(rho2/rho1) = {asym_model:.2f}"
          if hi > 0 else "Asymmetry: layer lies entirely below y = 0")
    print(f"Mean density in layer:  {rho_mean_layer:6.2f}   model sqrt(rho1 rho2) = {rho_model:.2f}"
          f"   (paper's simulation: 3.21)")
    print(f"Mean-flow KE removed:   {1 - ke / ke0:6.2f}   model 0.53, paper's simulation 0.58")
    print(f"Layer zero-momentum frame velocity: {v_frame:+.4f}")

    # Figure 6-style plot
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6), constrained_layout=True)
    curves = [
        (rho, rho0_i, r"$\langle\rho\rangle$"),
        (vx - v_frame, mom0_i / rho0_i - v_frame, r"$\langle v_x\rangle$ (layer frame)"),
        (0.5 * rho * (vx - v_frame) ** 2,
         0.5 * rho0_i * (mom0_i / rho0_i - v_frame) ** 2,
         r"$\frac{1}{2}\langle\rho\rangle\langle v_x\rangle^2$"),
    ]
    for ax, (now, init, label) in zip(axes, curves):
        ax.plot(y, init, color="0.6", ls="--", label="t = 0")
        ax.plot(y, now, color="k", label=tlabel)
        for edge in (lo, hi):
            ax.axvline(edge, color="tab:red", ls=":", lw=1)
        ax.set_xlim(lo - 0.3 * (hi - lo), hi + 0.3 * (hi - lo))
        ax.set_xlabel("y")
        ax.set_ylabel(label)
    axes[0].axhline(rho_model, color="tab:blue", lw=1, label=r"$\sqrt{\rho_1\rho_2}$")
    axes[0].legend(fontsize=8)
    fig.savefig("profiles.png", dpi=150)
    print("saved profiles.png")


if __name__ == "__main__":
    main()
