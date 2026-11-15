#!/bin/bash -l
# AGILE reduced 3D reference run of the Hillier & Arregui (2019) KH setup
# (80 x 208 x 400 cells, t = 0 -> 60) on one wICE GPU node (4 x P100).
# Needs in $VSC_DATA/kh-mixing/3d_ref_run:  mod_usr.fpp  Makefile  agile_ref3d.par
# Submit from the login node:   sbatch agile_ref3d.slurm
#SBATCH --export=NONE
#SBATCH --job-name=kh3d_ref
#SBATCH --output=%x_%j.out
#SBATCH --account=lp_edu_amrvac
#SBATCH --clusters=wice
#SBATCH --partition=gpu_p100
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=4
#SBATCH --time=04:00:00

# Modules from the course's AGILE example
module load Python/3.13
module load OpenMPI/5.0.8-nvidia-compilers-25.9-CUDA-12.9.1

export AGILE_DIR="$VSC_DATA/AGILE"
export PATH="$PATH:.:$AGILE_DIR:$AGILE_DIR/tools"
source "$AGILE_DIR/.venv/bin/activate"

echo "GPUs on this node:"
nvidia-smi -L

# Run in scratch: fast parallel output, room for ~16 GB of snapshots
SETUP="$VSC_DATA/kh-mixing/3d_ref_run"
RUNDIR="$VSC_SCRATCH/kh3d_ref"
mkdir -p "$RUNDIR/output"
cp "$SETUP"/{Makefile,mod_usr.fpp,agile_ref3d.par} "$RUNDIR"/
cd "$RUNDIR"

# Compile for the GPU (OpenACC), with the run's par file as build config
make -j $SLURM_CPUS_PER_TASK arch=nvidia OPENACC=1 CONFIG=agile_ref3d.par

# One MPI process per GPU
mpirun -np $SLURM_NTASKS ./agile -i agile_ref3d.par

echo "Done. Output in $RUNDIR/output"
