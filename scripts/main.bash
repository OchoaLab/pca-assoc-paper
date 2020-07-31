#################
### MAIN SIMS ###
#################

# all have: -k 10 -f 0.1
# this is like my older FST work
# but unlike my main gas-rgls sim (instead has: -k 3 -f 0.3)

# the large sample size simulation
Rscript pcs_test_LMM_PCA_n_1000_1.R # -n 1000 --m_causal 100 -g 1

# the small sample size simulation
Rscript pcs_test_LMM_PCA_n_1000_1.R -n 100 --m_causal 10 # -g 1

# the family structure simulation
Rscript pcs_test_LMM_PCA_n_1000_1.R -g 20 # -n 1000 --m_causal 100

# TODO:
# - separate steps more to run more parallel on cluster, save partial results
# - set up in a way that code can be reused for real datasets

#####################
### Human Origins ###
#####################

# make links to data
# based on: gas-rgls/scripts/data_links.bash
# use LD pruned only

# shared by R script calls further below
name=HoPacAll_ld_prune_1000kb_0.3

DATA_DIR='/home/viiia/dbs/humanOrigins'
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
# 1m25.639s viiiaX6

# get PCs in R, using my formula
# hope for similar average performance, although these PCs are very different than GCTA (except for top 2)
# NOTE: uses ROM version (known bias, does not match popular versions; I think this is best)
time Rscript real-01-pca-test.R --bfile $name
# 6m57.502s viiiaX6

# this creates auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
time Rscript real-02-subset-eigenvec.R --bfile $name
# 0m3.437s ideapad
# same but with standard PCA estimates (from my R code, kinship_std ROM version)
time Rscript real-02-subset-eigenvec.R --bfile $name --std

# calculates kinship matrix with popkin, to get mean kinship to pass to simtrait
time Rscript real-03-popkin.R --bfile $name
# 4m54.676s ideapad

# draws a random trait
time Rscript real-04-simtrait.R --bfile $name -r 1
# 0m2.166s ideapad

# loop for work computer
for rep in {2..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep
done

# GCTA runs, eventually do with all PCs
# runtime is remarkably constant!
time Rscript real-05-gcta.R --bfile $name -r 1 --n_pcs 0
# 5m36.827s ideapad
time Rscript real-05-gcta.R --bfile $name -r 1 --n_pcs 10
# 5m32.453s ideapad
time Rscript real-05-gcta.R --bfile $name -r 1 --n_pcs 90
# 5m37.847s ideapad

# loop that does all PCs in a given rep (local runs)
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs
    done
done
# same but with PCs run backwards (GCTA is only using half of processors on labby!, hopefully they meet in the middle)
for pcs in {90..0}; do
    for rep in {50..1}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs
    done
done


# PCA runs
time Rscript real-06-pca.R --bfile $name -r 1 --n_pcs 0
# 1m13.805s ideapad
time Rscript real-06-pca.R --bfile $name -r 1 --n_pcs 10
# 4m56.686s ideapad
time Rscript real-06-pca.R --bfile $name -r 1 --n_pcs 90
# 53m24.639s ideapad

# PCA runs (with plink)
time Rscript real-06-pca-plink.R --bfile $name -r 1 --n_pcs 0
# 0m4.355s ideapad
time Rscript real-06-pca-plink.R --bfile $name -r 1 --n_pcs 10
# 0m13.500s ideapad
time Rscript real-06-pca-plink.R --bfile $name -r 1 --n_pcs 90
# 3m45.629s ideapad

# loop that does all PCs in a given rep (local runs)
#for pcs in {0..90}; do
for pcs in {67..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs
    done
done
# same but with PCs run backwards (to run on another machine, hopefully they meet in the middle)
for pcs in {90..0}; do
    for rep in {50..1}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs
    done
done


# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90
# 38+8+5+2+1+3+5+2+1+1+1m ideapad (ran in parts, progressively)

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
# 0m32.963s ideapad (partial run)

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name
# 0m1.809s ideapad

# tests that p-value vectors have the right lengths of m_loci
# to make sure nothing was corrupted due to scripts stopping unexpectedly or incomplete file transfers
# (now test is peformed within real-07-auc-rmsd.R above too, but this can retest everything without recalculating expensive summary statistics).
time Rscript real-10-validate-pvals.R --bfile $name -r 50 --n_pcs 90
# 3m12.595s ideapad

###

# removes redundant, auxiliary GCTA PCA files
time Rscript real-02-subset-eigenvec.R --bfile $name --clean
# 0m0.473s ideapad
time Rscript real-02-subset-eigenvec.R --bfile $name --clean --std

###################
### SIMS REBOOT ###
###################

# let's rewrite the main simulations code to follow the real dataset analysis more in parallel, so that we take advantage of its embarrasingly parallel structure (for potential cluster runs)
# old code only ran R PCA anyway, for speed we need plink PCA (which requires plink BED files, etc), so lots of restructuring is needed regardless

### construct admixture proportions, ancestral inbreeding, and family structure if needed
# done once, shared across replicates

# the large sample size simulation
time Rscript sim-00-sim-pop.R # -n 1000 -g 1
# 0m1.047s ideapad

# the small sample size simulation
time Rscript sim-00-sim-pop.R -n 100 # -g 1
# 0m0.783s ideapad

# the family structure simulation
time Rscript sim-00-sim-pop.R -g 20 # -n 1000
# 0m13.656s ideapad (concurrent with plink, so probably would have been faster)

### construct genotype matrices

# the large sample size simulation
time Rscript sim-01-draw-geno.R -r 1 # -n 1000 -g 1
# 0m19.112s ideapad (concurrent with plink)

# hacks to use "real" data scripts on simulations
name="sim-n1000-k10-f0.1-s0.5-g1"

# draws a random trait
time Rscript sim-02-sim-trait.R --bfile $name -r 1
# 0m1.914s ideapad (concurrent with plink)

# preprocess with GCTA (makes GRM and max PCs)
time Rscript real-00-preprocess-gcta.R --bfile $name/rep-1
# 0m56.426s ideapad (concurrent with plink)

# get PCs in R, using my formula
# NOTE: uses ROM version (known bias, does not match popular versions; I think this is best)
# top 4 PCs match perfectly here, 5 agrees partially, after that it's a mess
time Rscript real-01-pca-test.R --bfile $name/rep-1
# 0m17.695s ideapad (concurrent with plink)

# this creates auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
time Rscript real-02-subset-eigenvec.R --bfile $name/rep-1
# 0m3.215s ideapad (concurrent with plink)
# same but with standard PCA estimates (from my R code, kinship_std ROM version)
time Rscript real-02-subset-eigenvec.R --bfile $name/rep-1 --std
# 0m3.422s ideapad (concurrent with plink)

# NOTE: unlike real dataset, here we don't need popkin estimates, we use true p_anc instead (for simulating trait earlier)!

# GCTA runs, eventually do with all PCs
# runtime is remarkably constant!
time Rscript real-05-gcta.R --sim --bfile $name -r 1 --n_pcs 0
# 0m40.506s ideapad (concurrent with plink)
time Rscript real-05-gcta.R --sim --bfile $name -r 1 --n_pcs 10
# 0m36.261s ideapad (concurrent with plink)
time Rscript real-05-gcta.R --sim --bfile $name -r 1 --n_pcs 90
# 0m37.381s ideapad (concurrent with plink)

# PCA runs (with plink)
time Rscript real-06-pca-plink.R --sim --bfile $name -r 1 --n_pcs 0
# 0m3.138s ideapad (concurrent with plink)
time Rscript real-06-pca-plink.R --sim --bfile $name -r 1 --n_pcs 10
# 0m3.094s ideapad (concurrent with plink)
time Rscript real-06-pca-plink.R --sim --bfile $name -r 1 --n_pcs 90
# 0m13.181s ideapad (concurrent with plink)

# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name

# tests that p-value vectors have the right lengths of m_loci
# to make sure nothing was corrupted due to scripts stopping unexpectedly or incomplete file transfers
# (now test is peformed within real-07-auc-rmsd.R above too, but this can retest everything without recalculating expensive summary statistics).
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90

###############
### JUNK??? ###
###############

# currently unused
# rmsd_auc_read: function to read rmsd or auc for pca or gcta (example attached)

# this series doesn't look fully functional
# pcs_test_n_100_trait: code to generate random trait and save all data
# pcs_test_n_100_pca: code to conduct association test for pca for one pc (pc need to be set in advance)
# pcs_test_n_100_gcta: code to conduct association test for gcta for one pc (pc need to be set in advance & this code still report the error that missing mala in both my desktop and laptop)
