# this script gathers available AUC and RMSD estimates (in separate, tiny files) into a big table

library(optparse)
library(readr)
library(dplyr)

# constants
methods <- c('pca-plink', 'gcta')
# output file name (big table)
file_table <- 'sum.txt'

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
                help = "Max replicates", metavar = "int")
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

# big table of interest
# initialize this way, it'll grow correctly
tib_main <- NULL

for ( rep in 1 : rep_max ) {
    # move higher to the "reps" location
    # this is so GCTA's temporary files don't overwrite files from other parallel runs
    dir_out <- paste0( 'rep-', rep )
    # skip reps that we haven't calculated at all
    if ( !dir.exists( dir_out ) )
        next
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
                    col_types = 'ciiddd'
                )
                # concatenate into bigger table
                tib_main <- bind_rows( tib_main, tib )
            }
        }
    }
    
    # move back down when done with this rep
    setwd( '..' )
}

# write the big table to file!
write_tsv(
    tib_main,
    file_table
)
