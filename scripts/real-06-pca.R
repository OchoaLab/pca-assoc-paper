# this scripts runs PCA on simulation, with a certain number of PCs

library(optparse)
library(readr)
library(BEDMatrix)
library(genio)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('gas_lm_optim.R')
source('gas_pca_optim.R')
setwd( dir_orig ) # go back to where we were

# the name is for dir only, actual file is just "data"
name_in <- 'data'
# the number of PCs in the file of precalculated values (just for retrieving that file)
n_pcs_max <- 90

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--n_pcs", type = "integer", default = 0,
                help = "Number of PCs to use", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 1,
                help = "Replicate number", metavar = "int")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep <- opt$rep
n_pcs <- opt$n_pcs

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# move higher to the "reps" location
# this is so GCTA's temporary files don't overwrite files from other parallel runs
dir_out <- paste0( 'rep-', rep )
setwd( dir_out )

###########
### PCA ###
###########

# genotypes, PCs are all in lower level (shared across reps)
name_in_lower <- paste0( '../', name_in )

# load with BEDMatrix
X <- BEDMatrix(
    name_in_lower,
    simple_names = TRUE
)
# load trait (in current level)
phen <- read_phen( name_in )
trait <- phen$pheno # extract column of interest

message("PCA")

if ( n_pcs == 0 ) {
    # run basic linear model without covariates
    obj <- gas_lm_optim(X, trait)
} else {
    # get PCs, incorporate as covariates

    # first load PCs file (all 90 of them)
    # these are our estimates (from kinship_std ROM version)
    # (also in lower level directory)
    pcs <- read_table2(
        paste0( name_in_lower, '-std-n_pcs_', n_pcs_max, '.eigenvec' ),
        col_names = FALSE
    )
    # remove first two columns (FAM and ID)
    pcs <- pcs[ , -(1:2) ]
    # now subset to keep only the number we wish to incorporate
    pcs <- pcs[ , 1 : n_pcs ]
    # pass as numeric matrix
    pcs <- as.matrix( pcs )
    
    # actual run
    obj <- gas_pca_optim(X, trait, pcs)
}

# all we care to preserve here are the p-values
pvals <- obj$pvals
# save into a file, simple human-readable format
write_lines(
    pvals,
    paste0( 'pvals_pca_', n_pcs, '.txt.gz' )
)
