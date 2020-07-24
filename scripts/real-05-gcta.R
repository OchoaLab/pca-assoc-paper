# this scripts runs GCTA on simulation, with a certain number of PCs

library(optparse)
library(readr)

# load new functions from external scripts
dir_orig <- getwd()
setwd("../../scripts") # scripts of main GAS project
source('gas_lmm_gcta.R')
source('paths.R')
setwd( dir_orig ) # go back to where we were

# the name is for dir only, actual file is just "data"
name_in <- 'data'
# GCTA: here there's no timing, so let's be as fast as possible (use all threads)
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

############
### GCTA ###
############

# file to create
file_out <- paste0( 'pvals_gcta_', n_pcs, '.txt.gz' )

# do not redo run if output was already present!
if ( file.exists( file_out ) )
    stop( 'Output already exists, skipping: ', file_out )

# GCTA with PCA
message("GCTA with ", n_pcs, " PCs")

# output name for GCTA runs should have number of PCs, so concurrent runs don't overwrite each other
name_out <- paste0( 'gcta_', n_pcs )

# actual GWAS
obj <- gas_lmm_gcta(
    gcta_bin = gcta_bin,
    name = paste0( '../', name_in ), # genotypes, GRM, PCA are all in lower level (shared across reps)
    name_phen = name_in, # phenotype is in current level
    name_out = name_out, # write outputs into current level, add number of PCs
    threads = threads,
    n_pcs = n_pcs
)

# all we care to preserve here are the p-values
pvals <- obj$pvals
# save into a file, simple human-readable format
write_lines(
    pvals,
    file_out
)

# clean up when we're done with gcta
# deletes GAS table only (mlma and log)
# files are in current level (matches earlier `name_out` option)
delete_files_gcta( name_out )
