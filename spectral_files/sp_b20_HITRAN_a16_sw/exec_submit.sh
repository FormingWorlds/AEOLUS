#!/bin/bash
#SBATCH --job-name=s20sw
#SBATCH --output=sp_b20_HITRAN_a16_sw.out
#SBATCH --error=sp_b20_HITRAN_a16_sw.err
#SBATCH --time=168:00:00
#SBATCH -A snic2022-1-1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
date
export OMP_NUM_THREADS=16
module load netCDF/4.3.2-HDF5-1.8.12-nsc1-intel-2018.u1-bare
source /home/x_ryabo/socrates_2211/socrates_main/set_rad_env
python sp_gen_script.py
date
