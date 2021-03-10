# reports actual m_causal values used in every simulation
# (for paper, noticed some very odd issues)

library(readr)
library(tibble)
library(optparse)

#################
### CONSTANTS ###
#################

file_simtrait <- 'simtrait.RData'
rep_max <- 50
# names of datasets
datasets <- c(
    'sim-n1000-k10-f0.1-s0.5-g1',
    'sim-n100-k10-f0.1-s0.5-g1',
    'sim-n1000-k10-f0.1-s0.5-g20',
    'HoPacAll_ld_prune_1000kb_0.3',
    #'hgdp_wgs_autosomes_ld_prune_1000kb_0.3',
    'hgdp_wgs_autosomes_ld_prune_1000kb_0.3_maf-0.01',
    #'all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1'
    'all_phase3_filt-minimal_ld_prune_1000kb_0.3_maf-0.01'
)
# directory name, needed in one mode (weird var name just stuck)
dir_phen <- 'const_herit_loci/'

############
### ARGV ###
############

# define options
option_list = list(
    make_option("--const_herit_loci", action = "store_true", default = FALSE, 
                help = "Causal coefficients constructed to result in constant per-locus heritability (saved in diff path)")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
const_herit_loci <- opt$const_herit_loci


############
### DATA ###
############

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

        # move in an additional level in this case (must exist)
        if ( const_herit_loci )
            setwd( dir_phen )
        
        # load file that has true m_causal
        load( file_simtrait )
        # and extract such value, save in vector
        m_causals[ rep ] <- length(causal_indexes)
        
        # go back down
        setwd( '..' )
        # move back an additional level in this case
        if ( const_herit_loci )
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
