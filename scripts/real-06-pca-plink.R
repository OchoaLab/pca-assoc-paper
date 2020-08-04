# this scripts runs PCA on simulation, with a certain number of PCs
# this version uses plink for the actual association, but uses my eigenvectors precomputed in R

library(optparse)
library(readr)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('gas_plink.R')
source('paths.R')
setwd( dir_orig ) # go back to where we were

# the name is for dir only, actual file is just "data"
name_in <- 'data'
# plink: here there's no timing, so let's be as fast as possible (use all threads)
threads <- 0

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
                help = "Replicate number", metavar = "int"),
    make_option("--sim", action = "store_true", default = FALSE, 
                help = "Genotypes are simulated (rather than real; alters location only)")
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

method <- 'pca-plink'

# message so we know where we're at
message(
    'rep: ', rep,
    ', method: ', method,
    ', pcs: ', n_pcs
)

# file to create
file_out <- paste0( 'pvals_', method, '_', n_pcs, '.txt.gz' )

# do not redo run if output was already present!
if ( file.exists( file_out ) )
    stop( 'Output already exists, skipping: ', file_out )

# genotypes, PCs:
# - in real data, are all in lower level (shared across reps)
# - in simulated data, are all in current level (not shared across reps)
name_in_lower <- if ( opt$sim ) name_in else paste0( '../', name_in )

# output name for plink runs should have number of PCs, so concurrent runs don't overwrite each other
name_out <- paste0( 'plink_', n_pcs )

file_covar <- paste0( name_in_lower, '-std-n_pcs_', n_pcs, '.eigenvec' )
# there's no covariates file to pass if we want zero PCs
# gas_plink will handle this correctly
if ( n_pcs == 0 )
    file_covar <- NULL

# actual run
obj <- gas_plink(
    plink2_bin,
    name = name_in_lower,
    name_phen = name_in,
    name_out = name_out, # write outputs into current level, add number of PCs
    file_covar = file_covar,
    threads = threads
)

# all we care to preserve here are the p-values
pvals <- obj$pvals
# save into a file, simple human-readable format
write_lines(
    pvals,
    file_out
)

# clean up when we're done with plink
# deletes GAS table and log only
# files are in current level (matches earlier `name_out` option)
delete_files_plink_gas( name_out )
