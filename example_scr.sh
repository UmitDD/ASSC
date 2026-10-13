#!/bin/bash
#SBATCH -J SS
#SBATCH --nodes=1
#SBATCH -t 2-00:00:00
#SBATCH -p compute

ulimit -s unlimited
module purge
export OMP_NUM_THREADS="1"

##########vasp-tasks##############
module load vasp/vasp.5.4.4.pl2-oneapi
export VASP="vasp_ncl"