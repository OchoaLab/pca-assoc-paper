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

###############
### JUNK??? ###
###############

# currently unused
# rmsd_auc_read: function to read rmsd or auc for pca or gcta (example attached)

# this series doesn't look fully functional
# pcs_test_n_100_trait: code to generate random trait and save all data
# pcs_test_n_100_pca: code to conduct association test for pca for one pc (pc need to be set in advance)
# pcs_test_n_100_gcta: code to conduct association test for gcta for one pc (pc need to be set in advance & this code still report the error that missing mala in both my desktop and laptop)
