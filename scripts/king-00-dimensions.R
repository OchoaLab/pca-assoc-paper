# version for king-cutoff data, no need for lots of detail here, let's keep it simple
# - checks FES only
# - compares m_causal values but does not store them in output table (they are trivial anyway)

library(readr)
library(tibble)
library(genio)

#################
### CONSTANTS ###
#################

file_simtrait <- 'simtrait.RData'
rep_max <- 50
verbose <- FALSE # for troubleshooting

# directory name, needed in one mode (weird var name just stuck)
dir_phen <- 'fes/'
# suffix to actual dirs
suffix <- '_king-cutoff-4'
# shared in all cases
reps <- 1 : rep_max

############
### DATA ###
############

# move to where the data is
setwd( '../data/' )

# read in some starter info
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )
# read precalculated dimensions data to get original sample sizes
dims <- read_tsv( 'dimensions.txt', col_types = 'cciicid' )

# use real ones only here
datasets <- datasets[ datasets$type == 'Real', ]
dims <- dims[ dims$type == 'Real', ]

# add original sample sizes, after checking for dataset agreement
stopifnot( all( datasets$name_paper == dims$name_paper ) )
datasets$n_ind_orig <- dims$n_ind

# add suffixes to dir names
datasets$name_dir <- paste0( datasets$name_dir, suffix )
# toss some data we don't need here and don't want in output
datasets$col <- NULL
datasets$lty <- NULL
datasets$type <- NULL # not needed after the above filtering step

# initialize other columns to avoid errors
datasets$m_loci <- NA
datasets$n_ind <- NA

# load each dataset
for ( i in 1 : nrow( datasets ) ) {
    # copy down name
    name_dir <- datasets$name_dir[ i ]
    if ( verbose )
        message( name_dir ) 
    # enter dir
    setwd( name_dir )
    
    # get data dimensions
    # for real datasets, the data is in this directory
    datasets$n_ind[ i ] <- count_lines( 'data.fam', verbose = verbose )
    datasets$m_loci[ i ] <- count_lines( 'data.bim', verbose = verbose )
    
    # data to save
    # gather for FES traits only, they should all agree!
    m_causals_fes <- vector( 'integer', rep_max )

    # navigate reps
    for ( rep in reps ) {
        # go into this other dir
        setwd( paste0( 'rep-', rep ) )

        # now read FES trait
        # move in an additional level in this case (must exist)
        setwd( dir_phen )
        # load file that has true m_causal
        load( file_simtrait )
        # and extract such value, save in vector
        m_causals_fes[ rep ] <- length(causal_indexes)
        
        # go down two levels now
        setwd( '../..' )
    }
    
    # merge lists
    m_causals <- m_causals_fes
    # if these don't all agree, complain and stop (expect complete agreement!)
    # error message includes complete list
    if ( length( unique( m_causals ) ) != 1 )
        stop( 'm_causals did not all agree: ', datasets$name_paper[ i ], ': ', toString( m_causals ) )
    
    # go back down
    setwd( '..' )
}

# calculate percent removals
datasets$ind_rm_pc <- round( ( datasets$n_ind_orig - datasets$n_ind ) / datasets$n_ind_orig * 100, 1 )

# clean up things before saving table, so it's ready for the paper to just load
# delete this column only
datasets$name_dir <- NULL
# original sample sizes are also redundant (appear in earlier table)
datasets$n_ind_orig <- NULL
# would rename other columns, but since we need LaTeX codes (math and superscripts), better to do that on LaTeX's side instead of here

# write table to a file
write_tsv( datasets, file = paste0( 'dimensions', suffix, '.txt' ) )
