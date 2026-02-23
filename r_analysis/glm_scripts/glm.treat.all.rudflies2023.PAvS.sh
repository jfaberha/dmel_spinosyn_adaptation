#!/bin/bash
#SBATCH --partition=cas	# Partition/Queue to use
#SBATCH --time=5-00:00:00
#SBATCH --job-name="glm.treat.all.2023.PAvS"
#SBATCH --cpus-per-task=10
#SBATCH --mem=400GB
#SBATCH --output=glm.treat.all.2023.PAvS.out
##SBATCH --mail-type=ALL	# Email notification: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=jfaberha@wsu.edu	# Email address for notifications

module load r
cd /scratch/user/jfaberha/20251208_111954/admera/gp_analysis/rudflies_2023_redo/r
Rscript --vanilla glm.treat.all.rudflies2023.PAvS.r
echo "Completed job on node $HOSTNAME"
