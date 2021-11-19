#!/bin/bash
#SBATCH --job-name=real-sim-draw-geno-%a
#SBATCH --output=real-sim-draw-geno-%a.out
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=1
##SBATCH --mail-user=alejandro.ochoa@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.0.0

name=hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_sim
rep=$SLURM_ARRAY_TASK_ID
time Rscript fit-04-draw-geno.R --bfile $name -r $rep --maf_real

module unload R/4.0.0
