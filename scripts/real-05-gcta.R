# this scripts runs GCTA on simulation, with a certain number of PCs

library(optparse)
library(readr)
library(genbin) # binary wrappers

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
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep <- opt$rep
n_pcs <- opt$n_pcs
threads <- opt$threads
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit

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

# new level to this hierarchy
if ( m_causal_fac != 10 ) {
    dir_out <- paste0( 'm_causal_fac-', m_causal_fac )
    # now move in there
    setwd( dir_out )
}

# new level to this hierarchy
if ( herit != 0.8 ) {
    dir_out <- paste0( 'h', herit )
    # now move in there
    setwd( dir_out )
}

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
# data is even lower still in this mode
if ( !opt$sim && m_causal_fac != 10 )
    name_in_lower <- paste0( '../', name_in_lower )
if ( !opt$sim && herit != 0.8 )
    name_in_lower <- paste0( '../', name_in_lower )

# path to GCTA PCs
file_covar <- paste0( name_in_lower, '-n_pcs_', n_pcs, '.eigenvec' )
# there's no covariates file to pass if we want zero PCs
# gcta_mlma will handle this correctly
if ( n_pcs == 0 )
    file_covar <- NULL

# adjust paths if using fes model
dir_phen <- '' # current dir
# use subdir instead in this case
if ( fes )
    dir_phen <- 'fes/'

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
data <- gcta_mlma(
    name_in_lower, # genotypes, GRM, PCA are all in lower level (shared across reps)
    name_phen = name_phen,
    name_out = name_out, # write outputs into current level, add number of PCs
    file_covar = file_covar,
    threads = threads
)

# all we care to preserve here are the p-values
pvals <- data$p
# save into a file, simple human-readable format
# NOTE: if pvals is NULL (GCTA returns that sometimes, when there's no model convergence/ model is underdetermined) then write_lines writes an empty file, but later read_lines doesn't like this!
# instead let's write the word NULL, so it's one non-empty line
write_lines(
    if ( is.null( pvals ) ) 'NULL' else pvals,
    file_out
)

# clean up when we're done with gcta
# deletes mlma and log
# files are in current level (matches earlier `name_out` option)
delete_files_gcta_mlma( name_out )
delete_files_log( name_out )
