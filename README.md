# Kelvin–Helmholtz mixing (Hillier & Arregui 2019 reproduction)

Team project for Computational Methods for Astrophysical Applications, KU Leuven.
We reproduce the 2D (x–y) version of the KH mixing-layer simulation in
Hillier & Arregui (2019, ApJ 885:101), Section 3, with MPI-AMRVAC, and test how
the numerical choices change the result.

## Repository layout

```
paper2d/
  mod_usr.t      initial condition (the paper's setup: rho 1/10, dV = 0.1, beta = 0.3, B_z,
                 random v_y seed, passive tracer). Change -> recompile with make.
  amrvac.par     REFERENCE RUN settings (hlld + cada3, 160x400 effective, t = 0..60).
                 Copy this file to make a variation; no recompile needed.
scripts/
  plot_temperature.py   temperature maps like the paper's Figure 5 (+ shared data loader)
  vy_growth.py          max|v_y| and mixing-layer width vs time (is KH growing?)
  profiles.py           x-averaged profiles like Figure 6 + the paper's key numbers
results/
  run_log.md     what each run was and what it showed -- add an entry for every run
figures/         plots worth sharing (PNGs here are committed)
```

## One-time setup

Needs MPI-AMRVAC installed with `AMRVAC_DIR` set, gfortran, OpenMPI and Python 3.
On Windows use WSL and keep everything in your Linux home (`~/...`), not under `/mnt/c/`.

```bash
git clone https://github.com/<username>/kh-mixing.git
cd kh-mixing
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Optional, hides harmless OpenMPI network warnings on a laptop/WSL:
`echo 'export OMPI_MCA_btl=^tcp' >> ~/.bashrc`

## Compile and run the reference run

```bash
cd paper2d
$AMRVAC_DIR/setup.pl -d=2 -v=3      # 2D grid, 3 vector components (needed for B_z)
make -j 4
mkdir -p output
mpirun -np 4 ./amrvac -i amrvac.par
```

About 20-40 min on 4 laptop cores. Check the startup printout lists 9 variables
(rho, m1, m2, m3, e, b1, b2, b3, trc1) and `seed type : random`.

## Analyse a run

From the repo root, with `.venv` active (`output/kh_ref` is the run's `base_filename`):

```bash
python scripts/vy_growth.py        paper2d/output/kh_ref
python scripts/plot_temperature.py paper2d/output/kh_ref 0 10 20 30   # t = 0, 20, 40, 60
python scripts/profiles.py         paper2d/output/kh_ref 20 30        # average t = 40-60
```

## Making a variation (resolution, limiter, solver, AMR, cores)

1. Copy the reference settings: `cp paper2d/amrvac.par paper2d/amrvac_<name>.par`
2. Change **only** the setting you are testing, plus `base_filename`
   (e.g. `'output/kh_koren'`) so you don't overwrite other runs.
3. Run with `mpirun -np 4 ./amrvac -i amrvac_<name>.par` (no recompile).
4. Run the three scripts on it and add an entry to `results/run_log.md`.

Common settings to vary (all in the par file):

| Assignment item | Setting | Block |
|---|---|---|
| 5 resolution | `domain_nx1`, `domain_nx2` | `&meshlist` |
| 6 limiter | `limiter` | `&methodlist` |
| 7 solver | `flux_scheme`, `time_stepper` | `&methodlist` |
| 8 AMR vs uniform | `refine_max_level` (1 = uniform) | `&meshlist` |
| 9 scaling | number after `-np` | command line |

## Repo conventions

- Snapshots (`*.dat`), build files and scratch plots are git-ignored.
- Figures worth sharing go in `figures/`; PNGs there are committed.
- Results text (`tee results/<run>_*.txt`) and `run_log.md` are committed.
- `git pull` before you start working, `git push` after you commit.
