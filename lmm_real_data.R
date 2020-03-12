# NOTES:
# - the time will also be precomputed this time, so we should consider its different values on different machines
# - LFA doesn't work with BEDMatrix.  Omitted for now
setwd("/dscrhome/yy222/gas-rgls-master/scripts1")
library(optparse)    # for terminal options
library(genio)       # to write BED files for external software
library(BEDMatrix)   # to load real genotypes with low memory usage
library(popkin)      # to estimate kinship
library(popkinsuppl) # for kinship_std (supports BEDMatrix and missing data, unlike the crappy local implementation)
library(readr)       # to write kinship matrix
library(tibble)      # to store data
library(lfa)         # GWAS gcatest
source("gas_plots.R")
source('paths.R')
source('gas_lmm_gemma.R')
source('gas_lmm_emmax.R')
source('gas_lmm_gcta.R')
source('kinship_to_evd.R')
# load new functions from external scripts
source('gas_lm_optim.R')
source('gas_pca_optim.R')
source('gas_rgls.R')


setwd("/dscrhome/yy222/human-differentiation-manuscript-master/data")
# define options
option_list = list(
  make_option(c("-i", "--input"), type = 'character', default = "human_origins_and_pacific_public", 
              help = "Path for input genotypes (BED/BIM/FAM files, excluding extension)"),
  make_option("--debug", action = "store_true", default = FALSE, 
              help = "debug mode (GCTA and GEMMA are fully verbose)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name_in <- opt$input
if (is.na(name_in))
  stop('`-i/--input` is mandatory!')
debug <- opt$debug

# save kinship data
file_kinship_data <- paste0(name_in, '.RData')

# HACK: before BEDMatrix, copy *-orig.fam to *.fam so it loads without issue
file.copy(
  paste0(name_in, '-orig.fam'),
  paste0(name_in, '.fam')
)

# load and process genotypess
# do once per run (these things don't change)
# load genotypes
message('BEDMatrix')
X <- BEDMatrix(name_in)

# estimate kinship, needed for accurate setting of heritability
# for benchmarks, as this is part of RGLS, let's invert kinship too and time the whole thing

# write TPED (for ancient EMMAX code)
message('plink_bed_to_tped')
plink_bed_to_tped(plink1_bin, name_in) # must be plink1 for now
fam <- read_fam(name_in)
kinship_estimate <- popkin(X, fam$fam)
kinship_inv_estimate <- solve(kinship_estimate)

kinship_estimate_old <- kinship_std(X) # estimate kinship the old way
eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old ) # get all eigenvalues

# compute eigenvectors from popkin estimate (for PCA "oracle/popkin" versions)
message('kinship_to_evd')
eigenvectors_estimate <- kinship_to_evd(kinship_estimate)

# write popkin kinship matrix in correct format
# shared by at least GEMMA and EMMAX
file_popkin_kin <- paste0(name_in, '.popkin_kinship.txt')
message('writing: ', file_popkin_kin)
kinship_tibble <- as_tibble( 2 * kinship_estimate, .name_repair = 'minimal')
write_tsv(kinship_tibble, file_popkin_kin, col_names = FALSE)

# files to delete when all is done
file_fam <- paste0(name_in, '.fam')
# test that it's actually there!
save(kinship_estimate, kinship_inv_estimate, file_popkin_kin, eigenvectors_estimate, eigenvectors_estimate_old, file = file_kinship_data)
