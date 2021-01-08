# this script gathers available AUC and RMSD estimates (in separate, tiny files) into a big table

library(optparse)
library(readr)
library(dplyr)

# constants
methods <- c('pca-plink', 'pca-plink-pure', 'gcta')
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
    make_option("--const_herit_loci", action = "store_true", default = FALSE, 
                help = "Causal coefficients constructed to result in constant per-locus heritability (saved in diff path)"),
    make_option("--m_causal_fac", type = "double", default = 10,
                help = "Proportion of individuals to causal loci", metavar = "double")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
name <- opt$bfile
rep_max <- opt$rep
n_pcs_max <- opt$n_pcs
archived <- opt$archived
const_herit_loci <- opt$const_herit_loci
m_causal_fac <- opt$m_causal_fac

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

# new level to this hierarchy
if ( m_causal_fac != 10 ) {
    dir_out <- paste0( 'm_causal_fac-', m_causal_fac )
    # now move in there
    setwd( dir_out )
}

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
    
    # move in an additional level in this case
    if ( const_herit_loci ) {
        dir_phen <- 'const_herit_loci/'
        # directory can be missing, skip in that case
        if ( !dir.exists( dir_phen ) ) {
            # just move down from rep-* case
            setwd( '..' )
            # skip rest
            next
        }
        # else move in
        setwd( dir_phen )
    }
    
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
    # move back an additional level in this case
    if ( const_herit_loci )
        setwd( '..' )
}

if ( archived ) {
    # just switch back to ordinary data location
    setwd( dir_orig )
    setwd( name )
    # new level to this hierarchy
    if ( m_causal_fac != 10 ) {
        dir_out <- paste0( 'm_causal_fac-', m_causal_fac )
        # now move in there
        setwd( dir_out )
    }
}

# if const_herit_loci is true, create a new directory (if it doesn't already exist) and move there
if ( const_herit_loci ) {
    dir_out <- 'const_herit_loci'
    # create first time, if needed
    if ( !dir.exists( dir_out ) )
        dir.create( dir_out )
    # move there
    setwd( dir_out )
}

# write the big table to file!
write_tsv(
    tib_main,
    file_table
)
