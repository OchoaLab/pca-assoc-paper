#!/bin/bash
#SBATCH --job-name=real-sim-%a
#SBATCH --output=real-sim-%a.out
#SBATCH --mem=32G
#SBATCH --ntasks-per-node=1
##SBATCH --mail-user=alejandro.ochoa@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.0.0
module load Plink/2.00a3LM

name=hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_sim
rep=$SLURM_ARRAY_TASK_ID
time Rscript fit-04-draw-geno.R --bfile $name -r $rep --maf_real
time Rscript sim-02-sim-trait.R --bfile $name -r $rep
time Rscript sim-02-sim-trait.R --bfile $name -r $rep --fes
time Rscript real-00-preprocess-gcta.R --bfile $name/rep-$rep
time Rscript real-01-pcs-plink.R --bfile $name/rep-$rep
time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep

module unload R/4.0.0
module unload Plink/2.00a3LM
