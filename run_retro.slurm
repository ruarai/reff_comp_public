#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --partition=cascade,sapphire
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --array=1-925
#SBATCH --job-name=periodicgpcast
#SBATCH -o ./slurm_out/periodicgpcast_dist_%a.txt # STDOUT

cd /data/gpfs/projects/punim1541/periodic_gp_cast/
  
module load R/4.5.0

export R_LIBS_USER="/home/ruarai/R_libs"

Rscript --vanilla R/forecast_renewal_retro.R ${SLURM_ARRAY_TASK_ID}

# Log this job's resource usage stats
my-job-stats -a -n -s