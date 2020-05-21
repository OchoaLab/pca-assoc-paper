library(optparse)

setwd("/dscrhome/yy222/gas-rgls-master/scripts1")
#load trait data
load("pcs_test_n_100_trait.RData")

# load new functions from external scripts

source('gas_lm_optim.R')
source('gas_pca_optim.R')
source('gas_lmm_gcta.R')
source('paths.R')
source('gas_plots.R')
# define options
option_list = list(
  make_option(c("-n", "--n_pc"), type = "integer", default =2, 
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
  
  
  # GCTA with zero PCs
  message("GCTA")
  gcta_bin      <- 'C:/Users/yiqiy/Downloads/gcta_1.92.4beta_win/bin/gcta64'
  # make GRM (shared by all GCTA runs)
  obj <- gas_lmm_gcta_kin(gcta_bin, name_out, debug = debug)
  # actual GWAS run
  obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads, debug = debug)
  rmsd_gcta_n_100_pcs<- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
  auc_gcta_n_100_pcs<- pvals_to_pr_auc(obj$pvals, causal_indexes)
  setwd(paste0("/hpc/group/biostat/yy222/LMM_PCA_rep_n_", n_ind,"/rep1"))
  
  write.table(rmsd_gcta_n_100_pcs, file = paste0("auc_pca_n_", n_ind, "_pcs_0.txt") )
  write.table(auc_gcta_n_100_pcs, file = paste0("auc_pca_n_", n_ind, "_pcs_0.txt") )
  
}else{
  
  # GCTA with PCA
  message("GCTAp")
  gcta_bin      <- 'C:/Users/yiqiy/Downloads/gcta_1.93.1beta_win/bin/gcta64'
  # gets PCs, saves them on a file
  gas_lmm_gcta_pca(gcta_bin, name_out, n_pcs = n_pc, debug = debug)
  
  # actual GWAS
  obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads, debug = debug, n_pcs = n_pc)
 
  
  M_rmsd_gcta <- pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
  M_auc_gcta <- pvals_to_pr_auc(obj$pvals, causal_indexes)
  write.table(rmsd_gcta_n_100_pcs, file = paste0("rmsd_gcta_n_", n_ind, "_pcs_",n_pc,".txt") )
  write.table(auc_gcta_n_100_pcs, file = paste0("auc_gcta_n_", n_ind, "_pcs_",n_pc,".txt") )
  
  # clean up when we're done with gcta
  delete_files_gcta(name_out) # deletes GAS table only, not GRM (as desired)
  # these get deleted separately
  delete_files_gcta_pca(name_out, n_pcs = n_pc)
}


