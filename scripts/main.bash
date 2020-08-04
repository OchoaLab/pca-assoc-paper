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

# draws random traits
# 0m2.166s ideapad
for rep in {1..50}; do
    time Rscript real-04-simtrait.R --bfile $name -r $rep
done

# GCTA runs
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-05-gcta.R --bfile $name -r $rep --n_pcs $pcs
    done
done

# PCA runs (with plink)
for pcs in {0..90}; do
    for rep in {1..50}; do
	time Rscript real-06-pca-plink.R --bfile $name -r $rep --n_pcs $pcs
    done
done

# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --bfile $name -r 50 --n_pcs 90
# 38+8+5+2+1+3+5+2+1+1+1+1+1+1m ideapad (ran in parts, progressively)
# 61m47.537s + 7m48.435s viiiaX6 (PCA complete, GCTA not yet)

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
# 0m14.780s ideapad (partial run)
# 0m26.445s viiiaX6 (partial)

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

# params shared across reps
# use right set for each case

## LARGE SAMPLE SIZE
# sample size (normal n=1000, small n=100)
n=1000
# number of causal loci should scale with n
mc=$(( n / 10 ))
# to create family structure (g=20) or not (g=1)
g=1
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

## FAMILY STRUCTURE
# sample size (normal n=1000, small n=100)
n=1000
# number of causal loci should scale with n
mc=$(( n / 10 ))
# to create family structure (g=20) or not (g=1)
g=20
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

## SMALL SAMPLE SIZE
# sample size (normal n=1000, small n=100)
n=100
# number of causal loci should scale with n
mc=$(( n / 10 ))
# to create family structure (g=20) or not (g=1)
g=1
# hacks to use "real" data scripts on simulations
name="sim-n$n-k10-f0.1-s0.5-g$g"

for rep in {1..50}; do
    # the large sample size simulation
    time Rscript sim-01-draw-geno.R -r $rep -n $n -g $g
    # 0m19.112s ideapad (concurrent with plink)

    # draws a random trait
    # only place where --m_causal gets specified (not in file paths, etc)
    time Rscript sim-02-sim-trait.R --bfile $name -r $rep --m_causal $mc
    # 0m1.914s ideapad (concurrent with plink)

    # preprocess with GCTA (makes GRM and max PCs)
    time Rscript real-00-preprocess-gcta.R --bfile $name/rep-$rep
    # 0m56.426s ideapad (concurrent with plink)

    # get PCs in R, using my formula
    # NOTE: uses ROM version (known bias, does not match popular versions; I think this is best)
    # top 4 PCs match perfectly here, 5 agrees partially, after that it's a mess
    time Rscript real-01-pca-test.R --bfile $name/rep-$rep
    # 0m17.695s ideapad (concurrent with plink)

    # this creates auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep
    # 0m3.215s ideapad (concurrent with plink)
    # same but with standard PCA estimates (from my R code, kinship_std ROM version)
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --std
    # 0m3.422s ideapad (concurrent with plink)

    # NOTE: unlike real dataset, here we don't need popkin estimates, we use true p_anc instead (for simulating trait earlier)!

    # GCTA runs
    # do all PCs of this rep
    for pcs in {0..90}; do
	time Rscript real-05-gcta.R --sim --bfile $name -r $rep --n_pcs $pcs
    done

    # # PCA runs (with plink)
    # do all PCs of this rep
    for pcs in {0..90}; do
	time Rscript real-06-pca-plink.R --sim --bfile $name -r $rep --n_pcs $pcs
    done
    
    # cleanup
    # remove redundant files only
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean
    time Rscript real-02-subset-eigenvec.R --bfile $name/rep-$rep --clean --std
done

# summarizes p-values into AUC and RMSD for each method/rep/pc
time Rscript real-07-auc-rmsd.R --sim --bfile $name -r 50 --n_pcs 90
# 46m10.111s viiiaX6 (large)
# 38m9.266s viiiaX6 (family)
# 33m3.880s viiiaX6 (small)

# read all individual summary tables (tiny files), gather into a master table
time Rscript real-08-table.R --bfile $name -r 50 --n_pcs 90
# 0m32.340s viiiaX6 (large)
# 0m32.287s viiiaX6 (family)
# 0m32.228s viiiaX6 (small)

# creates final plot for paper!
time Rscript real-09-figs.R --bfile $name

# tests that p-value vectors have the right lengths of m_loci
# to make sure nothing was corrupted due to scripts stopping unexpectedly or incomplete file transfers
# (now test is peformed within real-07-auc-rmsd.R above too, but this can retest everything without recalculating expensive summary statistics).
# --final requires that all files exist!
time Rscript real-10-validate-pvals.R --sim --bfile $name -r 50 --n_pcs 90 --final
# 2m54.797s - 3m43.697s viiiaX6

###

# a comparison of RMSD and lambda across all datasets
time Rscript real-11-inflation-across-datasets.R
