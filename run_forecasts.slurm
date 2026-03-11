#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=02:00:00
#SBATCH --partition=cascade,sapphire
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6G
#SBATCH --array=1-50
#SBATCH --job-name=periodicgpcast
#SBATCH -o ./slurm_out/periodicgpcast_dist_%a.txt # STDOUT

# export round_id="2025-05-08"

cd /data/gpfs/projects/punim1541/periodic_gp_cast/
  
module load R/4.2.1

export R_LIBS_USER="/home/ruarai/R_libs"

Rscript --vanilla R/forecast_renewal.R ${SLURM_ARRAY_TASK_ID} $round_id

# Log this job's resource usage stats
my-job-stats -a -n -s