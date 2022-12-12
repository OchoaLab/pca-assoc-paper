# this script gathers available AUC and RMSD estimates (in separate, tiny files) into a big table

library(optparse)
library(readr)
library(dplyr)

# constants
methods <- c('pca-plink-pure', 'gcta', 'gcta-labs')
# output file name (big table)
file_table <- 'sum.txt'
# dir for archived data, created manually on viiiaX6 only
dir_archive <- '~/dbs2/PCA/'

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
    make_option(c("-a", "--archived"), action = "store_true", default = FALSE, 
                help = "Data is archived (changes search paths)"),
    make_option("--fes", action = "store_true", default = FALSE, 
                help = "Use FES instead of RC trait model"),
    make_option("--herit", type = "double", default = 0.8, 
                help = "heritability", metavar = "double"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double"),
    make_option("--env1", type = "double", default = NA,
                help = "Variance of 1st (coarsest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double"),
    make_option("--env2", type = "double", default = NA,
                help = "Variance of 2nd (finest) level of environment (non-genetic) effects (default NA is no env)", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs
archived <- opt$archived
fes <- opt$fes
m_causal_fac <- opt$m_causal_fac
herit <- opt$herit
env1 <- opt$env1
env2 <- opt$env2

# do this consistency check early
if ( !is.na( env1 ) && is.na( env2 ) )
    stop( 'If --env1 is specified, must also specify --env2!' )

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )

if ( archived ) {
    dir_orig <- getwd()
    setwd( dir_archive )
}
setwd( name )

# remember this location, using global path, for easily returning across multiple levels when done
dir_base <- getwd()

# big table of interest
# initialize this way, it'll grow correctly
tib_main <- NULL

for ( rep in 1 : rep_max ) {
    # specify location of files to process, as many levels as needed
    dir_out <- paste0( 'rep-', rep, '/' )
    if ( fes )
        dir_out <- paste0( dir_out, 'fes/' )
    if ( m_causal_fac != 10 )
        dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
    if ( herit != 0.8 )
        dir_out <- paste0( dir_out, 'h', herit, '/' )
    if ( !is.na( env1 ) )
        dir_out <- paste0( dir_out, 'env', env1, '-', env2, '/' )

    # skip reps that we haven't calculated at all
    if ( !dir.exists( dir_out ) )
        next
    # else move to that destination
    setwd( dir_out )
    
    # start a big loop
    for ( method in methods ) {
        for ( n_pcs in 0 : n_pcs_max ) {
            # file to read
            file_sum <- paste0( 'sum_', method, '_', n_pcs, '.txt.gz' )

            # if output is already there, don't do anything (don't recalculate)
            if ( file.exists( file_sum ) ) {
                # read table
                tib <- read_tsv(
                    file_sum,
                    col_types = 'ciiddddddddddd'
                )
                # concatenate into bigger table
                tib_main <- bind_rows( tib_main, tib )
            }
        }
    }
    
    # return to base when done with this rep
    setwd( dir_base )
}

if ( archived ) {
    # just switch back to ordinary data location
    setwd( dir_orig )
    setwd( name )
}

# define destination directory
dir_out <- ''
if ( fes )
    dir_out <- paste0( dir_out, 'fes/' )
if ( m_causal_fac != 10 )
    dir_out <- paste0( dir_out, 'm_causal_fac-', m_causal_fac, '/' )
if ( herit != 0.8 )
    dir_out <- paste0( dir_out, 'h', herit, '/' )
if ( !is.na( env1 ) )
    dir_out <- paste0( dir_out, 'env', env1, '-', env2, '/' )

if ( dir_out != '' ) {
    # create first time, if needed
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out, recursive = TRUE )
    # move there
    setwd( dir_out )
}

# write the big table to file!
write_tsv(
    tib_main,
    file_table
)
