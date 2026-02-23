#!/bin/bash
#SBATCH --partition=cas	# Partition/Queue to use
#SBATCH --time=7-00:00:00
#SBATCH --job-name="glm.2023.PAvSvSEvE_leave1out"
#SBATCH --cpus-per-task=10
#SBATCH --mem=400GB
#SBATCH --output=glm.2023.PAvSvSEvE.leave1out.out
##SBATCH --mail-type=ALL	# Email notification: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=jfaberha@wsu.edu	# Email address for notifications

module load r
cd /scratch/user/jfaberha/20251201_112936/admera/gp_analysis/rudflies_2023_redo/r
Rscript --vanilla glm.rudflies2023.PAvSvSEvE.leave1out.r
echo "Completed job on node $HOSTNAME"
