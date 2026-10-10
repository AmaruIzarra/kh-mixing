# Kelvin–Helmholtz mixing (Hillier & Arregui 2019 reproduction)

Team project for Computational Methods for Astrophysical Applications, KU Leuven.

## One-time setup

Needs: MPI-AMRVAC installed with `AMRVAC_DIR` set, gfortran, OpenMPI, Python 3.
On Windows, use WSL and keep everything in your Linux home (`~/...`), not under `/mnt/c/`.

```bash
cd ~/codes
git clone https://github.com/AmaruIzarra/kh-mixing.git
cd kh-mixing
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

If `python3 -m venv` fails on Ubuntu: `sudo apt install python3.X-venv`
(the error message tells you the exact package name).

Optional, hides harmless OpenMPI network warnings on a laptop/WSL:
```bash
echo 'export OMPI_MCA_btl=^tcp' >> ~/.bashrc
```

## Compile and run the 2D setup

```bash
cd kh2d
$AMRVAC_DIR/setup.pl -d=2 -v=2
make -j 4
mpirun -np 4 ./amrvac -i amrvac.par
```

## Plot a snapshot

From the repo root, with `.venv` active:
```bash
python plot_snapshot.py kh2d/kh_2d0015.dat
```

## Repo conventions

- Snapshots (`*.dat`), build files and scratch plots are git-ignored.
- Figures worth sharing go in `figures/`. PNGs there are committed.
- `git pull` before you start working, `git push` after you commit.
