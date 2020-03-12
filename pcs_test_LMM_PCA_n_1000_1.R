setwd("/dscrhome/yy222/gas-rgls-master/scripts1")
library(optparse)

############
### ARGV ###
############

# define options
option_list = list(
  make_option(c("-p", "--plot"), action = "store_true", default = FALSE, 
              help = "create plot of existing data"),
  make_option("--simple", action = "store_true", default = FALSE, 
              help = "simple plot version (for a grant proposal)"),
  make_option(c("-l", "--lite"), action = "store_true", default = FALSE, 
              help = "run lite version (focuses on best and least redundant methods)"),
  make_option(c("-n", "--n_ind"), type = "integer", default = 1000, 
              help = "number of individuals", metavar = "int"),
  make_option(c("-m", "--m_loci"), type = "integer", default = 100000, 
              help = "number of loci", metavar = "int"),
  make_option(c("-k", "--k_subpops"), type = "integer", default = 10, 
              help = "admixture intermediate subpopulations", metavar = "int"),
  make_option(c("-f", "--fst"), type = "double", default = 0.1, 
              help = "FST (fixation index)", metavar = "double"),
  make_option(c("--bias_coeff"), type = "double", default = 0.5, 
              help = "admixture bias coeff", metavar = "double"),
  make_option(c("-g", "--generations"), type = "integer", default = 1, 
              help = "number of generations, for realistic local kinship", metavar = "int"),
  make_option("--herit", type = "double", default = 0.8, 
              help = "heritability", metavar = "double"),
  make_option("--m_causal", type = "integer", default = 100, 
              help = "num causal loci", metavar = "int"),
  make_option(c("-t", "--threads"), type = "integer", default = 1, 
              help = "number of threads (affects GCTA only)", metavar = "int"),
  make_option("--debug", action = "store_true", default = FALSE, 
              help = "debug mode (GCTA and GEMMA are fully verbose)"),
    make_option("--rep", type = "integer", default = 2, 
              help = "number of replication", metavar = "int")

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
n_ind <- opt$n_ind
m_loci <- opt$m_loci
k_subpops <- opt$k_subpops
fst <- opt$fst
bias_coeff <- opt$bias_coeff
generations <- opt$generations
m_causal <- opt$m_causal
herit <- opt$herit
lite <- opt$lite
threads <- opt$threads
simple <- opt$simple
debug <- opt$debug
rep<-opt$rep

# output path for BED files and figure
name_out <- paste0(
  'gas',
  '-n', n_ind,
  '-m', m_loci,
  '-k', k_subpops,
  '-f', fst,
  '-s', bias_coeff,
  '-mc', m_causal,
  '-h', herit,
  '-t', threads,
  '-g', generations
)
# mark simplified plots too



   

  
  # else run simulation
library(readr)       # to write kinship matrix
  library(tibble)      # to store data
  library(genio)       # to write BED files for external software
  library(popkin)      # to estimate kinship in RGLS
  library(popkinsuppl) # for PCA's kinship estimator
  library(lfa)         # GWAS gcatest
  library(gcatest)     # GWAS gcatest
  #library(qvalue)      # multiple hypothesis tests
  
  # standard code for a complex trait and an admixed population
  source('sim_geno_trait_k3.R')
  
  # load new functions from external scripts
  source('kinship_to_evd.R')
  source('gas_lm_optim.R')
  source('gas_pca_optim.R')
  source('gas_rgls.R')
  source('gas_lmm_gemma.R')
  source('gas_lmm_emmax.R')
  source('gas_lmm_gcta.R')
  source('gas_plots.R')
    style <- read_tsv('style.txt', col_types = 'cci')
    
  #Generate the matrix to store the value of
  M_rmsd_LM_n_100<-rep(NA,rep)
  M_auc_LM_n_100<-rep(NA,rep)
  M_rmsd_LMM_n_100<-matrix(0,rep,90)
  M_auc_LMM_n_100<-matrix(0,rep,90)
  M_rmsd<-matrix(0,rep,90)
  M_auc<-matrix(0,rep,90)
    M_auc_GCTA<-rep(NA,rep)
    M_rmsd_GCTA<-rep(NA,rep)

  ############
  ### SIMS ###
  ############
  for (i in 1:rep){
  obj <- sim_geno_trait_k3(
    n_ind = n_ind,
    m_loci = m_loci,
    m_causal = m_causal,
    k_subpops = k_subpops,
    bias_coeff = bias_coeff,
    generations = generations,
    herit = herit,
    verbose = TRUE,
    fst = fst
  )
  X <- obj$X
  admix_proportions <- obj$admix_proportions
  kinship <- obj$kinship
  trait <- obj$trait
  causal_indexes <- obj$causal_indexes
  
  ## how we'd get these values if we weren't setting them on command line
  ## n_ind <- ncol(X)
  ## m_loci <- nrow(X)
  ## k_subpops <- ncol(admix_proportions)
  ## sanity checks instead
  stopifnot( ncol(X) == n_ind )
  stopifnot( nrow(X) == m_loci )
  stopifnot( ncol(admix_proportions) == k_subpops )
  
  # write BED version for external code
  plink_data <- write_plink(name_out, X, pheno = trait, verbose = T)
  
  # write phenotype file
  # shared by EMMAX and GCTA
  write_phen(name_out, plink_data$fam, verbose = T)
  
    gcta_bin<-"/hpc/group/biostat/yy222/gcta_1.92.4beta1/gcta64"
  ########################
  ### ASSOCIATION TEST ###
  ########################

    # store runtimes in these tibbles
    pvals <- tibble(.rows = m_loci)
    
        #test that no need to repeat
         
       name <- "PCA"
        message(name)
       

  
      kinship_estimate_old <- kinship_std(X) # estimate kinship the old way
      eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old ) # get all eigenvalues

        
         name <- "GCTA"
    message(name)
    
    obj <- gas_lmm_gcta_kin(gcta_bin, name_out)
    #    file_gcta_kin <- obj$file

    
    obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads)
    pvals[[name]] <- obj$pvals
       M_rmsd_GCTA[i]<-pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
    M_auc_GCTA[i]<-pvals_to_pr_auc(obj$pvals, causal_indexes)

 name <- "LM"
    message(name)
    
    obj <- gas_lm_optim(X, trait)
   
    pvals[[name]] <- obj$pvals
     M_rmsd_LM_n_100[i]<-pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
    M_auc_LM_n_100[i]<-pvals_to_pr_auc(obj$pvals, causal_indexes)

        for (pcs in 1:90){
    ###########
    ### PCA ###
    ###########
    
    name <- "PCA"
    message(name)
       
    indexes <- 1 : pcs
    obj <- gas_pca_optim(X, trait, eigenvectors_estimate_old[, indexes ])
   
    pvals[[name]] <- obj$pvals
    M_auc[i,pcs]<-pvals_to_pr_auc(obj$pvals, causal_indexes)
    M_rmsd[i,pcs]<-pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd

    name <- "GCTA"
    message(name)
    
    obj <- gas_lmm_gcta_kin(gcta_bin, name_out)
    #    file_gcta_kin <- obj$file
  
    
    obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads)
    pvals[[name]] <- obj$pvals
   
    # GCTA with PCA

    name <- "GCTAp"
    message(name)
    
    # a new test, no timing this one
    # gets PCs, saves them on a file
    gas_lmm_gcta_pca(gcta_bin, name_out,  n_pcs = pcs,debug=debug)
    
     obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads, debug = debug, n_pcs = pcs)
    # clean up when we're done with gcta
    delete_files_gcta(name_out)
    # these get delted separately
    delete_files_gcta_pca(name_out, n_pcs = pcs)
    
    M_rmsd_LMM_n_100[i,pcs]<-pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
    M_auc_LMM_n_100[i,pcs]<-pvals_to_pr_auc(obj$pvals, causal_indexes)
    
        
          }
  
  }
  
  setwd("/dscrhome/yy222/LMM_PCA_rep_n_1000/rep1")
  write.table(M_auc, file="auc_k_fixed_k_10_pcs_1_90_n_1000_3_11_repeat_10_1.txt")
  write.table(M_rmsd, file="rmsd_k_fixed_k_10_pcs_1_90_n_1000_3_11_repeat_10_1.txt")
  write.table(M_auc_LMM_n_100, file="auc_k_fixed_k_10_pcs_1_90_n_1000_3_11_repeat_10_GCTAPC_1.txt")
  write.table(M_rmsd_LMM_n_100, file="rmsd_k_fixed_k_10_pcs_1_90_n_1000_3_11_repeat_10_GCTAPC_1.txt")
  write.table(M_auc_LM_n_100, file="auc_k_fixed_k_10_pcs_1_90_n_1000_3_11_repeat_10_LM_1.txt")
  write.table(M_rmsd_LM_n_100, file="rmsd_k_fixed_k_10_pcs_1_90_n_1000_3_11_repeat_10_LM_1.txt")
  write.table(M_auc_GCTA, file="auc_k_fixed_k_10_pcs_1_90_n_1000_3_11_repeat_10_GCTA_1.txt")
  write.table(M_rmsd_GCTA, file="rmsd_k_fixed_k_10_pcs_1_90_n_1000_3_11_repeat_10_GCTA_1.txt")
  

