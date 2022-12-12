# this script moves big intermediate files out of the normal workspace (which gets synced across machines)
# goal is to eventually store it on an inactive drive (no backups, not really needed)

library(optparse)
library(genio)

# constants
methods <- c('pca-plink-pure', 'gcta') # 'pca-plink', 
# move p-values (biggest) and individual summaries (should already be condensed into the master summary table)
bases <- c('pvals', 'sum') 
# destination dir, created manually on viiiaX6 only
dir_dest <- '~/dbs2/PCA/'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--bfile", type = "character", default = NA, 
                help = "Base name for input plink files (.BED/BIM/FAM)", metavar = "character"),
    make_option("--n_pcs", type = "integer", default = 90,
                help = "Max number of PCs", metavar = "int"),
    make_option(c("-r", "--rep"), type = "integer", default = 50,
                help = "Max replicates", metavar = "int"),
    make_option(c("-s", "--single_pc"), action = "store_true", default = FALSE,
                help = "Process a single PC (the value of `--n_pcs`) rather than the default of the full range (zero to the value of `--n_pcs`)"),
    make_option(c("-m", "--method"), type = "character", default = NA, 
                help = "Method to process (pca, lmm, or lmm-labs; default all)", metavar = "character"),
    make_option(c("-t", "--test"), action = "store_true", default = FALSE, 
                help = "Test run (makes sure no files are missing, but doesn't actually move anything.  However, output directories are created.)"),
    make_option(c("-u", "--unarchive"), action = "store_true", default = FALSE, 
                help = "Move files back from archive to original location."),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--env1", type = "double", default = NA,
                help = "Variance of 1st (coarsest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option("--env2", type = "double", default = NA,
                help = "Variance of 2nd (finest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option(c('-l', "--labs"), action = "store_true", default = FALSE, 
                help = "Include LMM with labels data")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2
labs <- opt$labs

# do this consistency check early
if ( !is.na( env1 ) && is.na( env2 ) )
    stop( 'If --env1 is specified, must also specify --env2!' )

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# add a third method in this case
if ( labs )
    methods <- c( methods, 'gcta-labs' )

# decide on range of PCs to process
pcs_vec <- if ( opt$single_pc ) n_pcs_max else 0 : n_pcs_max

# narrow down methods too if requested
if ( !is.na( opt$method ) ) {
    if ( opt$method == 'pca' ) {
        methods <- methods[1] # PCA is first one in original list (but actual name is more complicated)
    } else if ( opt$method == 'lmm' ) {
        methods <- methods[2] # LMM is second one in original list (but actual name is more complicated)
    } else if ( opt$method == 'lmm-labs' ) {
        methods <- methods[3]
    } else
        stop( 'Invalid `method` (must be "pca", "lmm", or not specified for both): ', method )
}

# move to where the data is
setwd( '../data/' )
setwd( name )

# remember this location, using global path, for easily returning across multiple levels when done
dir_base <- getwd()

# add same `name` dir to destination
dir_dest <- paste0( dir_dest, name )

if ( opt$unarchive ) {
    # if we want to unarchive, simply reverse source and destinations
    dir_base <- dir_dest
    dir_dest <- getwd()
    setwd( dir_base )
}

for ( rep in 1 : rep_max ) {
    # specify location of files to process, as many levels as needed
    dir_in <- paste0( 'rep-', rep, '/' )
    if ( fes )
        dir_in <- paste0( dir_in, 'fes/' )
    if ( m_causal_fac != 10 )
        dir_in <- paste0( dir_in, 'm_causal_fac-', m_causal_fac, '/' )
    if ( herit != 0.8 )
        dir_in <- paste0( dir_in, 'h', herit, '/' )
    if ( !is.na( env1 ) )
        dir_in <- paste0( dir_in, 'env', env1, '-', env2, '/' )
    # stop if missing
    if ( !dir.exists( dir_in ) )
        stop( 'whole rep ', rep, ' MISSING!' )
    setwd( dir_in )

    # create the same dir in the destination, if needed
    dir_out <- paste0( dir_dest, '/', dir_in )
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out, recursive = TRUE )
    
    # start a big loop
    for ( method in methods ) {
        # only one method is not expected to have PCs (important for final validation)
        pcs_vec_method <- if ( method == 'gcta-labs' ) 0 else pcs_vec
        for ( n_pcs in pcs_vec_method ) {
            # files to move
            for ( base in bases ) {
                file_in <- paste0( base, '_', method, '_', n_pcs, '.txt.gz' )

                if ( !file.exists( file_in ) )
                    # missing files are always fatal for this step
                    stop( 'MISSING: rep-', rep, ', ', file_in )
                # else move to destination
                file_out <- paste0( dir_out, '/', file_in )
                if ( file.exists( file_out ) )
                    stop( 'Output exists: rep-', rep, ', ', file_out )
                if ( !opt$test ) 
                    file.rename( file_in, file_out )
            }
        }
    }
    
    # return to base when done with this rep
    setwd( dir_base )
}

warns <- warnings()
if (!is.null( warns ) )
    print( warns )
