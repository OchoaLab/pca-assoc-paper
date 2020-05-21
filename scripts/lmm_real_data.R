setwd("/dscrhome/yy222/gas-rgls-master/scripts1")
#setwd("../scripts")
library(optparse)    # for terminal options
library(genio)       # read_fam
library(BEDMatrix)   # to load real genotypes with low memory usage
library(popkin)      # to estimate kinship (for simtrait)
library(popkinsuppl) # for kinship_std (supports BEDMatrix and missing data, unlike the crappy local implementation)
source('paths.R')
source('gas_lmm_gcta.R')
source('kinship_to_evd.R')

gcta_bin<-"/hpc/group/biostat/yy222/gcta_1.92.4beta1/gcta64"

setwd("/dscrhome/yy222/human-differentiation-manuscript-master/data")
#setwd('../../human-differentiation-manuscript/data/')
# define options
option_list = list(
  make_option(c("-i", "--input"), type = 'character', default = "human_origins_and_pacific_public", 
              help = "Path for input genotypes (BED/BIM/FAM files, excluding extension)")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name_in <- opt$input

# save kinship data
file_kinship_data <- paste0(name_in, '.RData')

message("GCTA")
gas_lmm_gcta_kin(gcta_bin, name_in)

# load and process genotypes
# do once per run (these things don't change)
message('BEDMatrix')
X <- BEDMatrix(name_in)
fam <- read_fam(name_in)

message('popkin')
kinship_estimate <- popkin(X, fam$fam)

message('kinship_std')
kinship_estimate_old <- kinship_std(X) # estimate kinship the old way
message('kinship_to_evd')
eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old ) # get all eigenvalues
# save PCA data
message('writing: ', file_kinship_data)
save(kinship_estimate, eigenvectors_estimate_old, file = file_kinship_data)

