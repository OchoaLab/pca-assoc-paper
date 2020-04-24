library(optparse)

setwd("/dscrhome/yy222/gas-rgls-master/scripts1")
#load trait data
load("pcs_test_n_100_trait.RData")

# load new functions from external scripts

source('gas_lm_optim.R')
source('gas_pca_optim.R')
source('paths.R')
source('gas_plots.R')
# define options
option_list = list(
  make_option(c("-n", "--n_pc"), type = "integer", default =1, 
              help = "number of principal components", metavar = "int"),
  make_option(c("-n", "--n_ind"), type = "integer", default = 100, # CHANGE: or 1000
              help = "number of individuals", metavar = "int"))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_pc <- opt$n_pc
n_ind <- opt$n_ind
########################
### ASSOCIATION TEST ###
########################

if ( n_pc==0) {
     #skip doing everything here }



message("LM")

obj <- gas_lm_optim(X, trait)
# NOTE zero PCS is position 1
 rmsd_pca_n_100_pcs<- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
 auc_pca_n_100_pcs<- pvals_to_pr_auc(obj$pvals, causal_indexes)

 write.table(rmsd_pca_n_100_pcs, file = paste0("auc_pca_n_", n_ind, "_pcs_0.txt") )
 write.table(auc_pca_n_100_pcs, file = paste0("auc_pca_n_", n_ind, "_pcs_0.txt") )

}else{

  message('PCs: ', n_pc)
  
  message("PCA")
  indexes <- 1 :  n_pc
  obj <- gas_pca_optim(X, trait, eigenvectors_estimate_old[, indexes ])
  
  rmsd_pca_n_100_pcs<- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
  auc_pca_n_100_pcs <- pvals_to_pr_auc(obj$pvals, causal_indexes)
  
  write.table(rmsd_pca_n_100_pcs, file = paste0("rmsd_pca_n_", n_ind, "_pcs_",n_pc,".txt") )
  write.table(auc_pca_n_100_pcs, file = paste0("auc_pca_n_", n_ind, "_pcs_",n_pc,".txt") )

}