# Run log

One entry per run: settings that differ from the reference, cost, key numbers, conclusion.
Generate the numbers with the three scripts and save them with `tee results/<run>_*.txt`.

## Reference run: hlld + cada3, paper's random seed (`paper2d/amrvac.par`)

- Setup: 2.5D MHD (2D x-y grid, 3 vector components), Hillier & Arregui (2019) Sec. 3.1:
  rho 1 (y<0) / 10 (y>0), tanh width 0.003, dV = 0.1 in the zero-momentum frame
  (V1 = 0.0909, V2 = -0.0091), p = 1/gamma = 0.6, beta = 0.3 -> B_z = 2,
  v_y seeded with uniform random noise |v_y| <= 0.01 dV in every cell, passive tracer.
- Box x in [-0.4, 0.4] (periodic), y in [-1.5, 0.5] (reflecting walls).
- Grid: 80x200 base, 2 AMR levels -> 160x400 effective (paper's x-y resolution), 20x20 blocks.
- Scheme: fourstep + hlld + cada3, CFL 0.8. Run to t = 60, snapshot every 2.
- dt ~ 9e-4, set by the fast magnetosonic speed (~2.2 in the light gas, 20x the shear).

Results:
- Growth: exponential phase t ~ 2-8 at rate ~0.5 per time unit (linear theory with a
  layer ~0.05 thick allows sigma <~ 0.57); saturation from t ~ 12 at max|v_y| ~ 0.04-0.07.
- Mixing-layer width (tracer 1-99%): ~0.06 during linear growth, then ~linear growth
  ~0.01 per time unit (~0.1 dV) to 0.60 at t = 60; plateaus at t ~ 28-34, 48-50 = vortex merging.
- Picture: many small vortices at t = 20 (like paper Fig. 5b); broad band of mixed gas
  at t = 60; 2D merging still enlarges structures late.
- Profiles (t = 40-60 average), compared with the paper:

  | Quantity | This run | Model | Paper's 3D sim |
  |---|---|---|---|
  | Layer edges / width | -0.403 / +0.157, 0.56 | - | width ~0.5 |
  | Asymmetry \|y_low\|/y_up | 2.56 | 3.16 | - |
  | Mean density in layer | 3.55 | 3.16 | 3.21 |
  | Mean-flow KE removed | 0.38 | 0.53 | 0.58 |

  Mass check: 72% of the layer on the light side -> 0.72*1 + 0.28*10 = 3.52.
- Interpretation: qualitative agreement (shift to the light side, mean flow vanishes
  in the layer centre where <rho> ~ sqrt(10), layer width ~ paper's). Remaining gap most
  likely from 2D: no 3D vortex disruption, no field-line bending (B_z only adds pressure).
- Files: results/ref_growth.txt, results/ref_profiles_t40-60.txt, figures/ref_*.png

## Lessons from earlier test runs (kept locally)

- tvdlf + minmod: the KH instability never develops at this resolution. Numerical
  diffusion (scaling with the fast speed ~2.2) widens the layer past the stability
  threshold before one e-folding time. Use hlld (or test this deliberately for item 7).
- A seed of 10 sine modes including the box-sized mode produces one box-sized vortex
  and a layer ~2x too wide; the paper's random noise is the reference choice.
