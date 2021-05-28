# reports actual m_causal values used in every simulation
# (for paper, noticed some very odd issues)
# checks both "inv" and "rand" traits (no need for separate runs)

library(readr)
library(tibble)
library(genio)

#################
### CONSTANTS ###
#################

file_simtrait <- 'simtrait.RData'
file_bnpsd <- 'bnpsd.RData'
rep_max <- 50

# directory name, needed in one mode (weird var name just stuck)
dir_phen <- 'const_herit_loci/'

############
### DATA ###
############

# move to where the data is
setwd( '../data/' )

# read in some starter info
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )
# toss some data we don't need here and don't want in output
datasets$col <- NULL
datasets$lty <- NULL

# initialize other columns to avoid errors
datasets$m_loci <- NA
datasets$n_ind <- NA
datasets$K <- NA
datasets$m_causal <- NA

# shared in all cases
reps <- 1 : rep_max

# load each dataset
for ( i in 1 : nrow( datasets ) ) {
    # copy down name
    name_dir <- datasets$name_dir[ i ]
    # enter dir
    setwd( name_dir )
    # decide if it's simulation or not, which decides where some genotype data ought to be
    is_sim <- grepl( 'sim', name_dir )
    
    # get data dimensions
    # for real datasets, the data is in this directory (not inside reps)
    if ( !is_sim ) {
        #datasets$n_ind[ i ] <- count_lines( 'data.fam', verbose = FALSE )
        datasets$m_loci[ i ] <- count_lines( 'data.bim', verbose = FALSE )
        # for K, have to read the annotation tables
        # read full fam data (unfortunately annotations have additional pops, so need this to subset)
        fam <- read_fam( 'data', verbose = FALSE )
        datasets$n_ind[ i ] <- nrow( fam )
        # read annotations
        subpop_info <- read_tsv('pops-annot.txt', comment = '#')
        # map subpopulations using sub-subpopulations
        fam$superpop <- subpop_info$superpop[ match( fam$fam, subpop_info$pop ) ]
        # now report K, this one is a range
        datasets$K[ i ] <- paste0(
            length( unique( fam$superpop ) ), # number of subpopulations is first number
            '-',
            length( unique( fam$fam ) ) # number of sub-subpopulations is last number
        )
    } else {
        # initialize these vectors (each rep has a different genotype matrix for full simulations)
        n_ind <- vector( 'integer', rep_max )
        m_loci <- vector( 'integer', rep_max )
        # get pop structure dimension from "bnpsd" data file
        load( file_bnpsd )
        # loads: admix_proportions, and either inbr_subpops or tree_subpops
        if ( exists( "inbr_subpops" ) ) {
            datasets$K[ i ] <- length( inbr_subpops )
            rm( inbr_subpops ) # remove for next round
        } else if ( exists( "tree_subpops" ) ) {
            datasets$K[ i ] <- length( tree_subpops$tip.label )
            rm( tree_subpops ) # remove for next round
        } else
            stop( 'Simulated data missing both `inbr_subpops` and `tree_subpops`!' )
    }
    
    # data to save
    # gather for both rand and inv (const_herit_loci) traits, they should all agree!
    m_causals_rnd <- vector( 'integer', rep_max )
    m_causals_inv <- vector( 'integer', rep_max )

    # navigate reps
    for ( rep in reps ) {
        # go into this other dir
        setwd( paste0( 'rep-', rep ) )

        # get dimensions for simulated data
        if ( is_sim ) {
            n_ind[ rep ] <- count_lines( 'data.fam', verbose = FALSE )
            m_loci[ rep ] <- count_lines( 'data.bim', verbose = FALSE )
        }

        # read "rand" trait first
        # load file that has true m_causal
        load( file_simtrait )
        # and extract such value, save in vector
        m_causals_rnd[ rep ] <- length(causal_indexes)

        # now read "inv" trait
        # move in an additional level in this case (must exist)
        setwd( dir_phen )
        # load file that has true m_causal
        load( file_simtrait )
        # and extract such value, save in vector
        m_causals_inv[ rep ] <- length(causal_indexes)
        
        # go down two levels now
        setwd( '../..' )
    }

    # check coherence of dimensions for simulations (where there's more than one rep)
    if ( is_sim ) {
        if ( length( unique( n_ind ) ) != 1 )
            stop( 'n_ind did not all agree: ', datasets$name_paper[ i ], ': ', toString( n_ind ) )
        if ( length( unique( m_loci ) ) != 1 )
            stop( 'm_loci did not all agree: ', datasets$name_paper[ i ], ': ', toString( m_loci ) )
        # store result in same tibble as other data
        # since all are now verified to be equal, just store first one
        datasets$n_ind[ i ] <- n_ind[ 1 ]
        datasets$m_loci[ i ] <- m_loci[ 1 ]
    }

    # merge lists
    m_causals <- c( m_causals_rnd, m_causals_inv )
    # if these don't all agree, complain and stop (expect complete agreement!)
    # error message includes complete list
    if ( length( unique( m_causals ) ) != 1 )
        stop( 'm_causals did not all agree: ', datasets$name_paper[ i ], ': ', toString( m_causals ) )
    # store result in same tibble as other data
    # since all are now verified to be equal, just store first one
    datasets$m_causal[ i ] <- m_causals[ 1 ]
    
    # go back down
    setwd( '..' )
}

# clean up things before saving table, so it's ready for the paper to just load
# delete this column only
datasets$name_dir <- NULL
# would rename other columns, but since we need LaTeX codes (math and superscripts), better to do that on LaTeX's side instead of here

# write table to a file
write_tsv( datasets, file = 'dimensions.txt' )
