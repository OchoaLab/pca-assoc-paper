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
# cleanup
rm ../data/$name/data.log
rm ../data/$name/data-n_pcs_90.log

# get PCs in R, using my formula
# hope for similar average performance, although these PCs are very different than GCTA (except for top 2)
# NOTE: uses ROM version (known bias, does not match popular versions; I think this is best)
time Rscript real-01-pca-test.R --bfile $name
# 6m57.502s viiiaX6

# this creates auxiliary GCTA PCA files (redundant, will delete when done with this analysis)
time Rscript real-02-subset-eigenvec.R --bfile $name
# 0m3.437s ideapad

# removes redundant, auxiliary GCTA PCA files
time Rscript real-02-subset-eigenvec.R --bfile $name --clean
# 0m0.473s ideapad

# calculates kinship matrix with popkin, to get mean kinship to pass to simtrait
time Rscript real-03-popkin.R --bfile $name
# 4m54.676s ideapad

###############
### JUNK??? ###
###############

# currently unused
# rmsd_auc_read: function to read rmsd or auc for pca or gcta (example attached)

# this series doesn't look fully functional
# pcs_test_n_100_trait: code to generate random trait and save all data
# pcs_test_n_100_pca: code to conduct association test for pca for one pc (pc need to be set in advance)
# pcs_test_n_100_gcta: code to conduct association test for gcta for one pc (pc need to be set in advance & this code still report the error that missing mala in both my desktop and laptop)