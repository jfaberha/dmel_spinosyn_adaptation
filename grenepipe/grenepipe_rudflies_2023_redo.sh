#!/bin/bash
#SBATCH --partition=cas # Partition/Queue to use
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=2
#SBATCH --job-name="gp_rudflies23redo"
#SBATCH --cpus-per-task=40
#SBATCH --mem=300GB
#SBATCH --output=gp_rudflies23redo.out
##SBATCH --mail-type=ALL        # Email notification: BEGIN,END,FAIL,ALL
##SBATCH --mail-user=jfaberha@wsu.edu   # Email address for notifications

module load miniconda3/3.12
source activate grenepipe_v0.13.2
cd /data/lab/rudman/gp_analysis/grenepipe-master
snakemake --use-conda --directory /scratch/user/jfaberha/20250811_113557/admera/gp_analysis/rudflies_2023_redo



