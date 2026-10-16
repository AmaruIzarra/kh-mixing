"""Temperature snapshots in the style of Hillier & Arregui (2019), Figure 5.

Usage (from the repo root, with .venv active):
    python scripts/plot_temperature.py paper2d/output/kh_paper 0 10 20 30

The first argument is the snapshot base name, the rest are snapshot numbers.
With dtsave_dat = 2, numbers 0, 10, 20, 30 are t = 0, 20, 40, 60.

Temperature is in the paper's non-dimensional form, T = gamma * p / rho,
so the light gas starts at T = 1 and the dense gas at T = 0.1.
"""
import sys

import matplotlib.pyplot as plt
import numpy as np
import yt

GAMMA = 5.0 / 3.0
YLIM = (-0.6, 0.4)  # same y-window as the paper's Figure 5

yt.set_log_level("warning")


def load_uniform(fname):
    """Load a snapshot and resample it onto a uniform grid at the finest level.

    Returns x, y (1D cell-centre coordinates), a dict of 2D arrays
    indexed [x, y], and the simulation time.
    """
    ds = yt.load(fname)
    level = ds.index.max_level
    dims = ds.domain_dimensions * ds.refine_by**level
    dims[2] = 1  # 2D run: a single cell in z
    cg = ds.covering_grid(level=level, left_edge=ds.domain_left_edge, dims=dims)

    # raw conserved variables, in code units
    f = {name: np.asarray(cg[("amrvac", name)])[:, :, 0]
         for name in ("rho", "m1", "m2", "m3", "e", "b1", "b2", "b3", "trc1")}

    rho = f["rho"]
    kinetic = 0.5 * (f["m1"]**2 + f["m2"]**2 + f["m3"]**2) / rho
    magnetic = 0.5 * (f["b1"]**2 + f["b2"]**2 + f["b3"]**2)
    p = (GAMMA - 1.0) * (f["e"] - kinetic - magnetic)

    f["p"] = p
    f["T"] = GAMMA * p / rho
    f["vx"] = f["m1"] / rho
    # Tracer as a fluid fraction (0 = light gas, 1 = dense gas). Depending on
    # the AMRVAC version it is stored either directly or as rho * tracer;
    # values well above 1 mean the latter.
    trc = f["trc1"]
    f["tracer"] = trc / rho if trc.max() > 1.5 else trc

    left = np.asarray(ds.domain_left_edge)
    right = np.asarray(ds.domain_right_edge)
    nx, ny = dims[0], dims[1]
    dx = (right[0] - left[0]) / nx
    dy = (right[1] - left[1]) / ny
    x = left[0] + dx * (np.arange(nx) + 0.5)
    y = left[1] + dy * (np.arange(ny) + 0.5)
    return x, y, f, float(ds.current_time)


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    base = sys.argv[1]
    numbers = [int(n) for n in sys.argv[2:]]

    fig, axes = plt.subplots(1, len(numbers), figsize=(3.2 * len(numbers), 4.2),
                             sharey=True, constrained_layout=True)
    axes = np.atleast_1d(axes)

    for ax, n in zip(axes, numbers):
        x, y, f, t = load_uniform(f"{base}{n:04d}.dat")
        im = ax.pcolormesh(x, y, f["T"].T, cmap="afmhot",
                           vmin=0.05, vmax=1.05, shading="auto")
        ax.set_ylim(*YLIM)
        ax.set_aspect("equal")
        ax.set_title(f"t = {t:.0f}")
        ax.set_xlabel("x")
    axes[0].set_ylabel("y")
    fig.colorbar(im, ax=axes, label=r"$T = \gamma p / \rho$", shrink=0.9)

    out = "temperature_panels.png"
    fig.savefig(out, dpi=150)
    print(f"saved {out}")


if __name__ == "__main__":
    main()
