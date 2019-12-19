setwd("/dscrhome/yy222/gas-rgls-master/scripts")
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
  make_option(c("-k", "--k_subpops"), type = "integer", default = 3, 
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
              help = "debug mode (GCTA and GEMMA are fully verbose)")
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



  style <- read_tsv('style.txt', col_types = 'cci')
  

  
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
  gcta_bin<-"'C:/Users/yiqiy/Documents/gcta_1.92.4beta2_win/bin/gcta64'"
  
  ############
  ### SIMS ###
  ############
  for (i in 1:2){
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
  
  #Generate the matrix to store the value of
  M_rmsd<-matrix(0,10,100)
  M_auc<-matrix(0,10,100)
  M_rmsd_LMM_n_1000_family<-matrix(0,10,100)
  M_auc_LMM_n_1000_family<-matrix(0,10,100)
  M_rmsd_LM_n_1000_family<-matrix(0,10,100)
  M_auc_LM_n_1000_family<-matrix(0,10,100)
  
  ########################
  ### ASSOCIATION TEST ###
  ########################

    # store runtimes in these tibbles
    times_kin <- tibble(.rows = 1)
    times_gas <- tibble(.rows = 1)
    pvals <- tibble(.rows = m_loci)
    
    for (pcs in 1:2){
    ###########
    ### PCA ###
    ###########
    
    name <- "PCA"
    message(name)
    
    times_kin[[name]] <- system.time({
      kinship_estimate_old <- kinship_std(X) # estimate kinship the old way
      eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old ) # get all eigenvalues
    })[3]
    
    indexes <- 1 : pcs
    obj <- gas_pca_optim(X, trait, eigenvectors_estimate_old[, indexes ])
    times_gas[[name]] <- obj$runtime
    pvals[[name]] <- obj$pvals
    
    name <- "GCTA"
    message(name)
    gcta_bin<-"C:/Users/yiqiy/Documents/gcta_1.92.4beta2_win/bin/gcta64"
    obj <- gas_lmm_gcta_kin(gcta_bin, name_out, debug = debug)
    #    file_gcta_kin <- obj$file
    times_kin[[name]] <- obj$runtime
    
    obj <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads, debug = debug)
    times_gas[[name]] <- obj$runtime
    pvals[[name]] <- obj$pvals
    
    # GCTA with PCA
    r <- pcs
    name <- "GCTAp"
    message(name)
    
    # a new test, no timing this one
    # gets PCs, saves them on a file
    gas_lmm_gcta_pca(gcta_bin, name_out, debug = debug, n_pcs = r)
    
    # don't record runtime here either
    pvals[[name]] <- gas_lmm_gcta(gcta_bin, name_out, m_loci = m_loci, threads = threads, debug = debug, n_pcs = r)$pvals
    
    # clean up when we're done with gcta
    delete_files_gcta(name_out)
    # these get delted separately
    delete_files_gcta_pca(name_out, n_pcs = r)
    
    
    name <- "LM"
    message(name)
    
    obj <- gas_lm_optim(X, trait)
    times_gas[[name]] <- obj$runtime
    pvals[[name]] <- obj$pvals
    
    #adjust the function to derive the rmsd and auc
    pvals_auc_collect  <- function (name_out, causal_indexes, pval_tab, times_kin, times_gas, style, cex_leg = 0.7) {
      # we shall generate other metrics in the plots below
      pval_rmsd <- list() # p-value mean distance from desired (under null)
      pr_auc <- list() # AUCs for P-R curves
      
      # reorder input tables as desired (style table)
      # also performs sanity checks (checks if style table is missing names)
      pval_tab <- reorder_tibble_names(style$name, pval_tab)
      times_kin <- reorder_tibble_names(style$name, times_kin)
      times_gas <- reorder_tibble_names(style$name, times_gas)
      
      # names and style order for p-value figs
      myNames <- names(pval_tab)
      order <- match(myNames, style$name)
      
      for (myName in myNames) {
        # sort p-values, compare to uniform quantiles, calculate RMSD
        obj <- pvals_to_null_rmsd(pval_tab[[myName]], causal_indexes)
        # store RMSD value
        pval_rmsd[[myName]] <- obj$rmsd
        # find row in style, to get col/lty from
        pr <- pvals_to_pr(pval_tab[[myName]], causal_indexes)
        # store AUCs
        pr_auc[[myName]] <- pr$auc.integral
        
        
      }
      return(cbind( pval_rmsd, pr_auc))
    }
    
    
    #store mean of p and auc values after removing the casual index
    p<-pvals_auc_collect(name_out, causal_indexes, pvals, times_kin, times_gas, style=style)
   
    
    M_auc[i,pcs]<-unlist(p[2,2])
    M_rmsd[i,pcs]<-unlist(p[2,1])
    M_rmsd_LMM_n_1000_family[i,pcs]<-unlist(p[4,2])
    M_auc_LMM_n_1000_family[i,pcs]<-unlist(p[4,2])
    M_rmsd_LM_n_1000_family[i,pcs]<-unlist(p[1,1])
    M_auc_LM_n_1000_family[i,pcs]<-unlist(p[1,2])
  }
  
  }

