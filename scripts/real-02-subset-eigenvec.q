#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=subset-eigenvec-%a
#SBATCH --output=subset-eigenvec-%a.out
##SBATCH --mem=1G
#SBATCH --ntasks-per-node=1
##SBATCH --mail-user=alejandro.ochoa@duke.edu
##SBATCH --mail-type=END,FAIL

module load R/4.0.0

#name=sim-n100-k10-f0.1-s0.5-g1
#name=sim-n1000-k10-f0.1-s0.5-g1
#name=sim-n1000-k10-f0.1-s0.5-g20
#name=hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01_geno-0.1_sim
#name=HoPacAll_ld_prune_1000kb_0.3_maf-0.01_sim
#name=tgp-nygc-autosomes_ld_prune_1000kb_0.3_maf-0.01_sim
rep=$SLURM_ARRAY_TASK_ID
time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --dcc --plink
time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --dcc

module unload R/4.0.0
