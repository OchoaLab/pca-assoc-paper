
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

