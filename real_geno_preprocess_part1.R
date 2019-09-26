

setwd("/dscrhome/yy222/gas-rgls-master/scripts")

library(optparse)    # for terminal options
library(genio)       # to write BED files for external software
library(BEDMatrix)   # to load real genotypes with low memory usage
library(popkin)      # to estimate kinship
library(popkinsuppl) # for kinship_std (supports BEDMatrix and missing data, unlike the crappy local implementation)
library(readr)       # to write kinship matrix
library(tibble)      # to store data
library(lfa)         # GWAS gcatest

source('kinship_to_evd.R')

# define options
option_list = list(
    make_option(c("-i", "--input"), type = 'character', default = NA, 
                help = "Path for input genotypes (BED/BIM/FAM files, excluding extension)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name_in <- "/dscrhome/yy222/gas-rgls-master/scripts/admix_g20"
if (is.na(name_in))
    stop('`-i/--input` is mandatory!')

# save kinship data
file_kinship_data <- paste0(name_in, '.RData')

# HACK: before BEDMatrix, copy *-orig.fam to *.fam so it loads without issue
file.copy(
    paste0(name_in, '-orig.fam'),
    paste0(name_in, '.fam')
)

# load and process genotypes
# do once per run (these things don't change)
# load genotypes
message('BEDMatrix')
X <- BEDMatrix(name_in)

# estimate kinship, needed for accurate setting of heritability
# for benchmarks, as this is part of RGLS, let's invert kinship too and time the whole thing

# store runtimes in these tibbles
times_kin <- tibble(.rows = 1)

name <- "RGLS"
message(name)

# in all cases, fam$fam are good labels
# load labels with genio
fam <- read_fam(name_in)

times_kin[[name]] <- system.time({
    kinship_estimate <- popkin(X, fam$fam)
    kinship_inv_estimate <- solve(kinship_estimate)
})[3]


name <- "PCA"
message(name)
times_kin[[name]] <- system.time({
    kinship_estimate_old <- kinship_std(X) # estimate kinship the old way
    eigenvectors_estimate_old <- kinship_to_evd( kinship_estimate_old ) # get all eigenvalues
})[3]


# compute eigenvectors from popkin estimate (for PCA "oracle/popkin" versions)
message('kinship_to_evd')
eigenvectors_estimate <- kinship_to_evd(kinship_estimate)

setwd("/dscrhome/yy222/real_data_set")

# save what we have to plot later!
save(kinship_estimate, kinship_inv_estimate, eigenvectors_estimate, eigenvectors_estimate_old, times_kin, file = file_kinship_data)
# , logistic_factors


