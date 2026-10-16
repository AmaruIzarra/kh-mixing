"""Track how fast the KH perturbation grows: max |v_y| and layer width per snapshot.

Usage (from the repo root, with .venv active):
    python scripts/vy_growth.py paper2d/output/kh_paper

Prints one line per snapshot and saves vy_growth.png (log scale).
If the instability is growing, max|v_y| should rise exponentially
(a straight line on the log plot) from ~1e-3 towards ~1e-2.
"""
import glob
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from plot_temperature import load_uniform  # noqa: E402


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    files = sorted(glob.glob(f"{sys.argv[1]}[0-9][0-9][0-9][0-9].dat"))
    if not files:
        sys.exit(f"no snapshots found matching {sys.argv[1]}NNNN.dat")

    times, vymax, widths = [], [], []
    print(f"{'t':>6} {'max|vy|':>10} {'layer width':>12}")
    for fname in files:
        x, y, f, t = load_uniform(fname)
        vy = f["m2"] / f["rho"]
        # layer width: y-range where the x-averaged tracer is between 1% and 99%
        tr = f["tracer"].mean(axis=0)
        mixed = y[(tr > 0.01) & (tr < 0.99)]
        width = mixed.max() - mixed.min() if mixed.size else 0.0

        times.append(t)
        vymax.append(np.abs(vy).max())
        widths.append(width)
        print(f"{t:6.1f} {vymax[-1]:10.2e} {width:12.3f}")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.5), constrained_layout=True)
    a1.semilogy(times, vymax, "o-")
    a1.set_xlabel("t")
    a1.set_ylabel(r"max $|v_y|$")
    a2.plot(times, widths, "o-")
    a2.set_xlabel("t")
    a2.set_ylabel("mixing-layer width (tracer 1-99%)")
    fig.savefig("vy_growth.png", dpi=150)
    print("saved vy_growth.png")


if __name__ == "__main__":
    main()
