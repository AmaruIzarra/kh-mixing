# Run log

## Run 1: reference, tvdlf + minmod (2026-10-16)
- Setup: paper2d (2.5D MHD, rho 1/10, dV = 0.1, beta = 0.3, B_z = 2)
- Grid: 80x200 base, 2 AMR levels (160x400 effective), blocks 20x20
- Scheme: fourstep + tvdlf + minmod, CFL 0.8
- Cost: 67,035 steps, 1243 s on 4 cores (1.38 core-hours); dt ~ 9e-4 (fast magnetosonic CFL)
- Timing: 50% ghost cells, 47% computing (small blocks -> communication-heavy)
- Result: NO KH vortices. max|vy| decays (3.4e-3 -> 2e-4).
  Layer width 0.005 -> 0.175 by t = 2 (numerical diffusion), 0.575 by t = 60.
- Interpretation: tvdlf diffusion scales with the fast speed (~2.2), 20x the shear (0.1).
  Layer exceeds kδ ~ 1 for the longest box mode (k = 2π/0.8, δ ~ 0.13) before one
  e-fold of growth (~4.4 time units), so no unstable mode fits in the box.
- Files: results/tvdlf_minmod_growth.txt, figures/tvdlf_minmod_*.png

## Run 2: hlld + cada3 (2026-10-16)
- Same setup and grid as Run 1; only flux_scheme (tvdlf -> hlld) and limiter (minmod -> cada3) changed
- Result: KH develops. max|vy| grows 3.4e-3 -> ~6e-2 by t ~ 12 (rate ~0.35 between t=2-4), saturates ~0.06-0.1 (~dV).
  Layer width 0.055 at t=2 (vs 0.175 tvdlf), ~0.96 by t = 60 (wider than box: 2D merging, geometric limit).
- Conclusion: solver choice decides whether the instability exists at all at this resolution.
- Files: results/hlld_cada3_growth.txt, figures/hlld_cada3_*.png
- Profiles (t = 40-60 average): edges y = -0.642 / +0.278 (width 0.92, ~2x paper's)
  asymmetry 2.32 (model 3.16), mean density 3.72 (model 3.16, paper 3.21),
  mean-flow KE removed 0.32 (model 0.53, paper 0.58).
  Mass check: 0.698*1 + 0.302*10 = 3.72, consistent with measured density.
  Matches paper qualitatively: shift to light side, mean KE -> 0 in layer centre (y ~ -0.13)
  where <rho> ~ sqrt(10). Quantitatively off: one box-sized 2D vortex dominates.
