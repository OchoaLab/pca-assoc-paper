# this script moves big intermediate files out of the normal workspace (which gets synced across machines)
# goal is to eventually store it on an inactive drive (no backups, not really needed)

library(optparse)
library(genio)

# constants
methods <- c('pca-plink', 'pca-plink-pure', 'gcta')
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
    make_option(c("-t", "--test"), action = "store_true", default = FALSE, 
                help = "Test run (makes sure no files are missing, but doesn't actually move anything.  However, output directories are created.)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs

# stop if name is missing
if ( is.na(name) )
    stop('`--bfile` terminal option is required!')

# move to where the data is
setwd( '../data/' )
setwd( name )

# create same `name` dir in output, if needed
dir_dest <- paste0( dir_dest, name )
if ( !dir.exists( dir_dest ) )
    dir.create( dir_dest )

for ( rep in 1 : rep_max ) {
    # move higher to the "reps" location
    dir_rep <- paste0( 'rep-', rep )
    setwd( dir_rep )

    # create the same dir in the destination, if needed
    dir_out <- paste0( dir_dest, '/', dir_rep )
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out )
    
    # start a big loop
    for ( method in methods ) {
        for ( n_pcs in 0 : n_pcs_max ) {
            # files to move
            for ( base in bases ) {
                file_in <- paste0( base, '_', method, '_', n_pcs, '.txt.gz' )

                if ( !file.exists( file_in ) )
                    # missing files are always fatal for this step
                    stop( 'MISSING: ', file_in )
                # else move to destination
                file_out <- paste0( dir_out, '/', file_in )
                if ( !opt$test ) 
                    file.rename( file_in, file_out )
            }
        }
    }
    
    # move back down when done with this rep
    setwd( '..' )
}
