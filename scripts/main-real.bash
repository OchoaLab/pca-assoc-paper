# make links to data
# based on: gas-rgls/scripts/data_links.bash
# use LD pruned only

# shared by R script calls further below
name='HoPacAll_ld_prune_1000kb_0.3'
DATA_DIR='/home/viiia/dbs/humanOrigins'

# version for HGDP WGS
#name='hgdp_wgs_autosomes_ld_prune_1000kb_0.3'
name='hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01'
DATA_DIR='/home/viiia/dbs/hgdp_wgs'

# version for TGP
#name='all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1'
name='all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01'
DATA_DIR='/home/viiia/dbs/tgp/plink2'

# shared steps (for each in turn)
cd ../data/
mkdir $name
cd $name

ln -s "$DATA_DIR/$name.bed" data.bed
ln -s "$DATA_DIR/$name.bim" data.bim
ln -s "$DATA_DIR/$name.fam" data.fam

# return to scripts dir
cd ../../scripts/
# preprocess with GCTA (makes GRM and max PCs)
time Rscript real-00-preprocess-gcta.R --bfile $name
# 1m25.639s viiiaX6 HO
# 1m26.897s ideapad HGDP
# 1m51.228s ideapad HGDP MAF
# 19m50.614s ideapad TGP
# 2m27.177s ideapad TGP thinned
# 2m3.016s ideapad TGP MAF

# get PCs in R, using my formula
# hope for similar average performance, although these PCs are very different than GCTA (except for top 2)
# NOTE: uses ROM version (known bias, does not match popular versions; I think this is best)
# time Rscript real-01-pca-test.R --bfile $name # Not repeated on TGP after m_causal bugs were fixed and thinning performed
# 6m57.502s viiiaX6 HO
# 5m21.411s ideapad HGDP
# 30m45.184s labbyDuke TGP # ideapad and viiiaX6 ran out of mem (~16G RAM)

# get PCs using plink2
time Rscript real-01-pcs-plink.R --bfile $name
# 0m39.102s ideapad HO
# 0m25.582s ideapad HGDP
# 0m18.244s ideapad HGDP MAF
# 2m48.309s ideapad TGP
# 0m5.659s ideapad TGP thinned, removed kinship_std comparison
# 0m32.804s ideapad TGP MAF

# calculates kinship matrix with popkin, to get mean kinship to pass to simtrait
time Rscript real-03-popkin.R --bfile $name
# 4m54.676s ideapad HO
# 8m37.050s ideapad HGDP
# 3m57.744s ideapad HGDP MAF
# 255m35.418s ideapad TGP
# 20m25.597s ideapad TGP thinned
# 21m53.466s ideapad TGP MAF

# draws random traits
# 0m2.166s ideapad HO (each rep)
# 0m21.066s ideapad TGP (each rep)
# 0m2.014s ideapad TGP thinned (each rep)
for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep
done

# create auxiliary PCA files from plink2 (redundant, will delete when done with this analysis)
#time Rscript real-02-subset-eigenvec.R --bfile $name --std
time Rscript real-02-subset-eigenvec.R --bfile $name --plink
# PCA runs (with pure plink)
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs --plink
    done
done
# removes redundant, auxiliary plink2 PCA files
#time Rscript real-02-subset-eigenvec.R --bfile $name --clean --std
time Rscript real-02-subset-eigenvec.R --bfile $name --clean --plink

# create auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
time Rscript real-02-subset-eigenvec.R --bfile $name
# 0m3.437s ideapad
# GCTA runs
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs
    done
done
# removes redundant, auxiliary GCTA PCA files
time Rscript real-02-subset-eigenvec.R --bfile $name --clean
# 0m0.473s ideapad

# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90
# all are plink + GCTA, labbyDuke 12 threads
# 19m16.850s HO
# 174m22.306s HGDP
# 71m22.370s TGP
# 80m32.881s TGP MAF

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
# 0m20.941s labbyDuke HO
# 0m27.567s labbyDuke HGDP
# 0m24.078s labbyDuke TGP
# 0m24.390s labbyDuke TGP MAF

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name
# 0m1.809s ideapad
# # OBSOLETE for partial runs, to plot only complete reps (best for slowest test: TGP)
# time Rscript real-09-figs.R --bfile $name --complete
# # OBSOLETE compares PCAs only (for internal purposes only); Not redone after m_causal bug was found
# time Rscript real-09-figs.R --bfile $name --pca

# a summary of "best PCs" in a simple analysis
Rscript real-13-stats.R --bfile $name
# Human Origins
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90    81 # OLD: 90,80
# 2 pca-plink-pure auc       34    33 # OLD: 34,31
# 3 gcta           rmsd       5     0
# 4 gcta           auc       12     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)
#
# HGDP MAF
#   method         metric  best   min
# 1 pca-plink-pure rmsd      57    26
# 2 pca-plink-pure auc       19     8
# 3 gcta           rmsd       0     0
# 4 gcta           auc        3     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)
#
# TGP MAF
#   method         metric  best   min
# 1 pca-plink-pure rmsd      68    39
# 2 pca-plink-pure auc       11     6
# 3 gcta           rmsd       0     0
# 4 gcta           auc       10     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)
#
# HGDP
#   method         metric  best   min
# 1 pca-plink-pure rmsd      35    17 # OLD: 72,17
# 2 pca-plink-pure auc       11     4
# 3 gcta           rmsd       0     0
# 4 gcta           auc        0     0
# best rmsd: pca-plink-pure (tie) # OLD: significant
# best auc: gcta (significant) # OLD tie
#
# TGP
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90     6 # OLD: 4,4
# 2 pca-plink-pure auc        7     4
# 3 gcta           rmsd       4     0
# 4 gcta           auc        5     0
# best rmsd: pca-plink-pure (significant)
# best auc: gcta (tie)

# tests that p-value vectors have the right lengths of m_loci
# to make sure nothing was corrupted due to scripts stopping unexpectedly or incomplete file transfers
# (now test is peformed within real-07-auc-rmsd.R above too, but this can retest everything without recalculating expensive summary statistics).
# --final requires that all files exist!
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90 --final
# 5m12.610s labbyDuke HO
# 35m17.690s labbyDuke HGDP
# 12m26.386s labbyDuke HGDP MAF
# 13m33.181s labbyDuke TGP
# 15m10.066s labbyDuke TGP MAF

# archive p-values and individual summary files (move out of space that gets synced between computers; these files take up loads of space and are no longer needed)
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90

# a comparison of RMSD and lambda across all datasets
time Rscript real-11-inflation-across-datasets.R
# model fit:
# rmsd ~ a * (lambda^b - 1) / (lambda^b + 1)
#         a         b 
# 0.5607461 0.6221887 
# log-linear approx: log(lambda) = RMSD * 5.73
# threshold map (sigmoidal): lambda = 1.05, RMSD = 0.00851
# threshold map (log-linear): lambda = 1.05, RMSD = 0.00851
#
# OLD
#         a         b 
# 0.5481480 0.6381526
# log-linear approx: log(lambda) = RMSD * 5.72
# threshold map (sigmoidal): lambda = 1.05, RMSD = 0.00853
# threshold map (log-linear): lambda = 1.05, RMSD = 0.00853

# reports on actual m_causal values used in all sims (used for paper, and to catch an unexpected error from an early run!)
time Rscript real-14-report-m-causal.R
# sim-n1000-k10-f0.1-s0.5-g1: 100
# sim-n100-k10-f0.1-s0.5-g1: 10
# sim-n1000-k10-f0.1-s0.5-g20: 100
# HoPacAll_ld_prune_1000kb_0.3: 292
# hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01: 93
# all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01: 250

# final plot gathers all three datasets into a single multipanel figure
time Rscript real-15-plots-big.R --real

### MAF code
## first gather MAF vectors from raw data
time Rscript real-16-mafs.R
# 1m57.285s ideapad
# 2m18.001s viiiaX6 MAF
## then make plot
time Rscript real-17-mafs-plot.R
# 0m13.937s ideapad

########################
### const_herit_loci ###
########################

# addition to complement existing analysis with alternative trait from const_herit_loci model
# not sure if it will make a meaningful difference
# NOTE: many steps that depend on genotypes only aren't redone (are shared from prev run)

for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep --const_herit_loci
done

# NOTE: these actually run on DCC through batch.R (requires manual edits)

# redo if needed
time Rscript real-02-subset-eigenvec.R --bfile $name --dcc --plink
# PCA runs (with pure plink)
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs --plink --const_herit_loci
    done
done
time Rscript real-02-subset-eigenvec.R --bfile $name --dcc --clean --plink

# redo if needed
time Rscript real-02-subset-eigenvec.R --bfile $name --dcc
# GCTA runs
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs --const_herit_loci
    done
done
time Rscript real-02-subset-eigenvec.R --bfile $name --dcc --clean

time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90 --const_herit_loci
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --const_herit_loci
time Rscript real-09-figs.R --bfile $name --const_herit_loci
Rscript real-13-stats.R --bfile $name --const_herit_loci
# Human Origins
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90    71 # OLD: 90,65
# 2 pca-plink-pure auc       34    16 # OLD: 34,9
# 3 gcta           rmsd       0     0
# 4 gcta           auc        1     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)
#
# HGDP MAF
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90    31
# 2 pca-plink-pure auc       17    15
# 3 gcta           rmsd       0     0
# 4 gcta           auc        1     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)
#
# TGP MAF
#   method         metric  best   min
# 1 pca-plink-pure rmsd      51    34
# 2 pca-plink-pure auc        9     8
# 3 gcta           rmsd       0     0
# 4 gcta           auc        1     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)
#
# HGDP
#   method         metric  best   min
# 1 pca-plink-pure rmsd      90    87 # OLD 2,2
# 2 pca-plink-pure auc        4     1
# 3 gcta           rmsd       0     0
# 4 gcta           auc        1     0
# best rmsd: pca-plink-pure (significant)
# best auc: gcta (significant) # OLD: tie
#
# TGP
#   method         metric  best   min
# 1 pca-plink-pure rmsd       1     0 # OLD 0,0
# 2 pca-plink-pure auc        2     2
# 3 gcta           rmsd       0     0
# 4 gcta           auc        0     0
# best rmsd: pca-plink-pure (significant)
# best auc: gcta (significant)

time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90 --final --const_herit_loci
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --const_herit_loci -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 --const_herit_loci

time Rscript real-15-plots-big.R --real --const_herit_loci

# a comparison of RMSD and lambda across all datasets
time Rscript real-11-inflation-across-datasets.R --const_herit_loci
# model fit:
# rmsd ~ a * (lambda^b - 1) / (lambda^b + 1)
#         a         b 
# 0.5291012 0.6772195 
# log-linear approx: log(lambda) = RMSD * 5.58
# threshold map (sigmoidal): lambda = 1.05, RMSD = 0.00874
# threshold map (log-linear): lambda = 1.05, RMSD = 0.00874

# reports on actual m_causal values used in all sims (used for paper, and to catch an unexpected error from an early run!)
time Rscript real-14-report-m-causal.R --const_herit_loci
# sim-n1000-k10-f0.1-s0.5-g1: 100
# sim-n100-k10-f0.1-s0.5-g1: 10
# sim-n1000-k10-f0.1-s0.5-g20: 100
# HoPacAll_ld_prune_1000kb_0.3: 292
# hgdp_wgs_autosomes_ld_prune_1000kb_0.3: 93
# all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1: 250


####################
### m_causal_fac ###
####################

# NOTE: this is a failed/abandoned experiment

#fac=100
#fac=70
#fac=60
#fac=55
fac=50
#fac=45
#fac=40

for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep --m_causal_fac $fac --const_herit_loci
done
# 0m4.383s ideapad HGDP each
# 0m3.740s ideapad TGP each

# PCA version
time Rscript real-02-subset-eigenvec.R --bfile $name --plink
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs --plink --m_causal_fac $fac --const_herit_loci
    done
done
# 1m1.205s ideapad HGDP
# 0m27.351s ideapad TGP
time Rscript real-02-subset-eigenvec.R --bfile $name --plink --clean

# GCTA version
time Rscript real-02-subset-eigenvec.R --bfile $name
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs --m_causal_fac $fac --const_herit_loci
    done
done
time Rscript real-02-subset-eigenvec.R --bfile $name --clean
# 8m51.743s ideapad HGDP
# 20m42.459s ideapad TGP

# ends up doing just the one file
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90 --m_causal_fac $fac --const_herit_loci
# 0m19.992s ideapad HGDP (both runs)
# 0m6.597s ideapad TGP (both runs)

time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90 --m_causal_fac $fac --const_herit_loci
time Rscript real-09-figs.R --bfile $name --m_causal_fac $fac --const_herit_loci

###################
### FIT FOR SIM ###
###################

# we already had popkin estimate because we needed params for simtrait (all real datasets)

# what follows has been run for TGP only!

# manually copied annotations file from Storey Lab project
# for TGP
cp ~/docs/ochoalab/storey/fst/simulations/all_phase3_filt-minimal/pops-annot.txt ../data/$name/
# for HGDP
cp ~/docs/ochoalab/storey/fst/simulations/hgdp_wgs_autosomes/pops-annot.txt ../data/$name/
# for HO
cp ~/docs/ochoalab/storey/fst/simulations/ho16ep/pops-annot.txt ../data/$name/

# this script calculates subpopulations kinship matrix, which we'll fit tree to
# also produces a nice plot that validates the estimate
time Rscript fit-01-popkin-subpops.R --bfile $name
# 4s ideapad
# NOTE: label/margin sizes tuned for TGP, terrible for HGDP and HO

# fit tree to real data, visualize success
time Rscript fit-02-tree.R --bfile $name

# set up simulation params Q and Psi tree
time Rscript fit-03-sim-pop.R --bfile $name

# from here on, the real data name needs a '_sim' suffix
# the rest of the pipeline is the same as for the regular simulations
name=$name"_sim"
for rep in {1..50}; do
    time Rscript fit-04-draw-geno.R --bfile $name -r $rep

    time Rscript sim-02-sim-trait.R --bfile $name -r $rep

    time Rscript real-00-preprocess-gcta.R --bfile $name/rep-$rep

    time Rscript real-01-pcs-plink.R --bfile $name/rep-$rep

    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean

    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs --plink
    done
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --plink --clean
done

time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90
# 14m22.146s # ideapad 4 threads
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
time Rscript real-09-figs.R --bfile $name
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --final
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90 -t # test first!
time Rscript real-12-archive-pvals.R --bfile $name -r 50 --n_pcs 90
Rscript real-13-stats.R --bfile $name
#   method         metric  best   min
# 1 pca-plink-pure rmsd      81    18
# 2 pca-plink-pure auc       15    14
# 3 gcta           rmsd       0     0
# 4 gcta           auc        1     0
# best rmsd: gcta (significant)
# best auc: gcta (significant)

# NOTE: this is not beta-p-inv!!!
