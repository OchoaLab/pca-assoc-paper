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
                help = "Genotypes are simulated (rather than real; alters paths only)"),
    make_option("--dcc", action = "store_true", default = FALSE, 
                help = "Duke Compute Cluster runs (alters paths only)"),
    make_option(c("-t", "--threads"), type = "integer", default = 0, 
                help = "number of threads (default use all cores if not DCC, 1 core if DCC)", metavar = "int"),
    make_option("--const_herit_loci", action = "store_true", default = FALSE, 
                help = "Causal coefficients constructed to result in constant per-locus heritability (saved in diff path)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep <- opt$rep
n_pcs <- opt$n_pcs
threads <- opt$threads
const_herit_loci <- opt$const_herit_loci

# figure out what threads should be
# above default should be modified if --dcc and if the threads aren't zero (default, means use all, which we should never do on a cluster!)
if ( opt$dcc && threads == 0 )
    threads <- 1

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
if ( opt$dcc ) {
    # on DCC we go here, no name for simplicity (data is only temporarily there, so I won't have more than one dataset there at the time)
    setwd( '/work/ao128/' )
} else {
    setwd( '../data/' )
}
setwd( name )

# move higher to the "reps" location
# this is so GCTA's temporary files don't overwrite files from other parallel runs
dir_out <- paste0( 'rep-', rep )
setwd( dir_out )

############
### GCTA ###
############

method <- 'gcta'

# message so we know where we're at
message(
    'rep: ', rep,
    ', method: ', method,
    ', pcs: ', n_pcs
)

# genotypes, PCs:
# - in real data, are all in lower level (shared across reps)
# - in simulated data, are all in current level (not shared across reps)
name_in_lower <- if ( opt$sim ) name_in else paste0( '../', name_in )

# adjust paths if using const_herit_loci model
dir_phen <- '' # current dir
# use subdir instead in this case
if ( const_herit_loci )
    dir_phen <- 'const_herit_loci/'

# only these are in dir_phen
name_phen <- paste0( dir_phen, name_in )
# output name for GCTA runs should have number of PCs, so concurrent runs don't overwrite each other
name_out <- paste0( dir_phen, 'gcta_', n_pcs )
# file to create
file_out <- paste0( dir_phen, 'pvals_', method, '_', n_pcs, '.txt.gz' )

# do not redo run if output was already present!
if ( file.exists( file_out ) )
    stop( 'Output already exists, skipping: ', file_out )

# actual GWAS
obj <- gas_lmm_gcta(
    gcta_bin = gcta_bin,
    name = name_in_lower, # genotypes, GRM, PCA are all in lower level (shared across reps)
    name_phen = name_phen,
    name_out = name_out, # write outputs into current level, add number of PCs
    threads = threads,
    n_pcs = n_pcs
)

# all we care to preserve here are the p-values
pvals <- obj$pvals
# save into a file, simple human-readable format
# NOTE: if pvals is NULL (GCTA returns that sometimes, when there's no model convergence/ model is underdetermined) then write_lines writes an empty file, but later read_lines doesn't like this!
# instead let's write the word NULL, so it's one non-empty line
write_lines(
    if ( is.null( pvals ) ) 'NULL' else pvals,
    file_out
)

# clean up when we're done with gcta
# deletes GAS table only (mlma and log)
# files are in current level (matches earlier `name_out` option)
delete_files_gcta( name_out )
