# reports actual m_causal values used in every simulation
# (for paper, noticed some very odd issues)

library(readr)
library(tibble)

# constants
file_simtrait <- 'simtrait.RData'
rep_max <- 50
# names of datasets
datasets <- c(
    'sim-n1000-k10-f0.1-s0.5-g1',
    'sim-n100-k10-f0.1-s0.5-g1',
    'sim-n1000-k10-f0.1-s0.5-g20',
    'HoPacAll_ld_prune_1000kb_0.3',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3',
    'all_phase3_filt-minimal_ld_prune_1000kb_0.3'
)

# move to where the data is
setwd( '../data/' )

# shared in all cases
reps <- 1 : rep_max

# load each dataset
for ( dataset in datasets ) {
    # enter dir
    setwd( dataset )

    # data to save
    m_causals <- vector( 'integer', rep_max )

    # navigate reps
    for ( rep in reps ) {
        # go into this other dir
        setwd( paste0( 'rep-', rep ) )

        # load file that has true m_causal
        load( file_simtrait )
        # and extract such value, save in vector
        m_causals[ rep ] <- length(causal_indexes)
        
        # go back down
        setwd( '..' )
    }

    # if these all agree, make a summarized report
    if ( length( unique( m_causals ) ) == 1 ) {
        message( dataset, ': ', unique( m_causals ) )
    } else {
        # this includes complete list
        message( dataset, ': ', toString( m_causals ) )
    }
    
    # go back down
    setwd( '..' )
}
