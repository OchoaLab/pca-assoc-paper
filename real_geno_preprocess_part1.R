# for terminal options
library(optparse)

# should be a command line option, meh
verbose <- TRUE

############
### ARGV ###
############
setwd("/dscrhome/yy222/gas-rgls-master/scripts1")
# define options
option_list = list(
  make_option(c("-i", "--input"), type = 'character', default = "human_origins_and_pacific_public", 
              help = "Path for input genotypes (BED/BIM/FAM files, excluding extension)"),
  make_option(c("-p", "--plot"), action = "store_true", default = FALSE, 
              help = "create plot of existing data"),
  make_option("--simple", action = "store_true", default = FALSE, 
              help = "simple plot version (for a grant proposal)"),
  make_option("--debug", action = "store_true", default = FALSE, 
              help = "debug mode (GCTA and GEMMA are fully verbose)"),
  make_option(c("-l", "--lite"), action = "store_true", default = FALSE, 
              help = "run lite version (focuses on best and least redundant methods)"),
  make_option("--herit", type = "double", default = 0.8, 
              help = "heritability", metavar = "double"),
  make_option("--m_causal", type = "integer", default = 100, 
              help = "num causal loci", metavar = "int"),
  make_option(c("-t", "--threads"), type = "integer", default = 1, 
              help = "number of threads (affects GCTA only)", metavar = "int")
)


 # load new functions from external scripts
    source('gas_lm_optim.R')
    source('gas_pca_optim.R')
    source('gas_rgls.R')
    source('gas_lmm_gemma.R')
    source('gas_lmm_emmax.R')
    source('gas_lmm_gcta.R')
    source('paths.R')
    source('gas_plots.R')
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name_in <- opt$input
if (is.na(name_in))
  stop('`-i/--input` is mandatory!')
m_causal <- opt$m_causal
herit <- opt$herit
lite <- opt$lite
threads <- opt$threads
simple <- opt$simple
debug <- opt$debug

# output path for BED files and figure
# NOTE: here name_in may (and does) have slashes because it's in a directory, so this only works if gas-* is not a prefix as in the original script
name_out <- paste0(
  name_in,
  '-gas',
  '-mc', m_causal,
  '-h', herit,
  '-t', threads
)
# mark simplified plots too
# NOTE: again use suffixes only, name_out has directory slashes
if (lite)
  name_out <- paste0(name_out, '-lite')

# saved data from earlier run (for easier plotting)
file_data <- paste0(name_out, '.RData')

  
  # else run simulation
  
  library(BEDMatrix) # to load real genotypes with low memory usage
  library(simtrait)  # to simulate complex trait
  library(readr)     # to write kinship matrix
  library(tibble)    # to store data
  library(genio)     # to write BED files for external software
  library(gcatest)   # GWAS gcatest
  
rep<-2
  
  ############
  ### SIMS ###
  ############
  setwd("/dscrhome/yy222/human-differentiation-manuscript-master/data")
  # load precomputed kinship matrix, etc:
  # kinship_estimate, kinship_inv_estimate, file_popkin_kin, eigenvectors_estimate, eigenvectors_estimate_old, times_kin, file_emmax_kin, file_gemma_kin
  # OMITTED FOR NOW: , logistic_factors
  file_kinship_data <- paste0(name_in, '.RData')
  if (verbose)
    message('loading: ', file_kinship_data)
  load(file_kinship_data)
  
  # HACK: before BEDMatrix, copy *-orig.fam to *.fam so it loads without issue
  file.copy(
    paste0(name_in, '-orig.fam'),
    paste0(name_in, '.fam')
  )
  
  # load genotypes
  if (verbose)
    message('BEDMatrix')
  X <- BEDMatrix(name_in)
  
  # in all cases, fam$fam are good labels
  # load data with genio
  fam <- read_fam(name_in, verbose = verbose)
  
  M_rmsd<-matrix(0,rep,90)
  M_auc<-matrix(0,rep,90)
  M_rmsd_gcta<-matrix(0,rep,90)
  M_auc_gcta<-matrix(0,rep,90)
  M_rmsd_LM<-rep(NA,rep)
  M_auc_LM<-rep(NA,rep)
  
  for (i in 1:2){
  # simulate trait!
  # can be done inside an inner loop to simulate thousands of traits per genotype dataset
  if (verbose)
    message('simtrait')
  # NOTE: using kinship version, only choice for real data
  obj_trait <- sim_trait(
    X = X,
    m_causal = m_causal,
    herit = herit,
    kinship = kinship_estimate
  )
  trait <- obj_trait$trait
  causal_indexes <- obj_trait$causal_indexes
  
  # data dimensions
  # simulation is only BEDMatrix, but meh
  m_loci <- if (class(X) == 'BEDMatrix') ncol(X) else nrow(X)
  
  # add phenotype to FAM table
  # safe since initial phenotypes were all 1 (essentially no data)
  fam$pheno <- trait
  # write phenotype file
  # shared by EMMAX and GCTA
  write_phen(name_in, fam, verbose = verbose)
  # for GEMMA (save to overwrite, original is safe elsewhere)
  write_fam(name_in, fam, verbose = verbose)
  
  
  
  # store runtimes in these tibbles
  times_gas <- tibble(.rows = 1)
  pvals <- tibble(.rows = m_loci)
  name <- "LM"
  message(name)
  
  obj <- gas_lm_optim(X, trait)
  
  pvals[[name]] <- obj$pvals
  
  M_rmsd_LM[i]<-pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
  M_auc_LM<-pvals_to_pr_auc(obj$pvals, causal_indexes)
  
  for (pcs in 1:2){
    
    
    name <- "PCA"
    message(name)
    
    indexes <- 1 : pcs
    obj <- gas_pca_optim(X, trait, eigenvectors_estimate_old[, indexes ])
    M_rmsd[i,pcs]<-pvals_to_null_rmsd(obj$pvals, causal_indexes)$rmsd
    M_auc[i,pcs]<-pvals_to_pr_auc(obj$pvals, causal_indexes)
    
  ############
  ### GCTA ###
  ############
  
  name <- "GCTA"
   gcta_bin<-"/hpc/group/biostat/yy222/gcta_1.92.4beta1/gcta64"
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
  }}

  
  write.table(M_auc_gcta, file="auc_gcta.txt")
  write.table(M_rmsd_gcta, file="rmsd_gcta.txt")
  write.table(M_auc, file="auc_pca.txt")
  write.table(M_rmsd, file="rmsd_pca.txt")
  write.table(M_auc_LM, file="auc_lm.txt")
  write.table(M_rmsd_LM, file="rmsd_lm.txt")
  