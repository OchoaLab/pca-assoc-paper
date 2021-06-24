# diagnosing issues by looking at trait normality and MAF of causal loci

library(readr)
library(ochoalabtools)
library(simtrait) # for rmsd()
library(BEDMatrix)

# constants
rep_max <- 50

#################
### FUNCTIONS ###
#################

trait_normality_rmsd <- function( trait ) {
    # standardize data so measurement to expected quantiles is easier
    trait <- drop( scale( trait ) )
    # sort too
    trait <- sort( trait )

    # now get expected quantiles for a standard normal!
    n <- length( trait )
    trait_exp <- qnorm( ppoints( n ) )

    # now compute a simple RMSD from these two, store immediately
    return( rmsd( trait, trait_exp ) )
}

causal_maf_min <- function( X, causal_indexes ) {
    # calculate MAFs of causal loci
    # use BEDMatrix orientation (loci_on_cols = TRUE)!
    mafs <- colMeans( X[ , causal_indexes ], na.rm = TRUE ) / 2
    # let's store only the minimum per rep
    return( min( mafs, na.rm = TRUE ) )
}

# processes both trait types together for greatest efficiency
# in addition to calculating trait normality, also gets MAF distributions of causal loci
process_reps <- function( is_sim ) {
    # data we are collecting
    rmsds_rand <- vector( 'numeric', rep_max )
    rmsds_inv <- vector( 'numeric', rep_max )
    mafs_rand <- vector( 'numeric', rep_max )
    mafs_inv <- vector( 'numeric', rep_max )
    
    if ( !is_sim )
        X <- BEDMatrix( 'data' )
    
    for ( rep in 1 : rep_max ) {
        # move to where the data is
        setwd( paste0( 'rep-', rep ) )
        
        # genotype data is here always (for sims), not in const_herit_loci
        if ( is_sim )
            X <- BEDMatrix( 'data' )

        # process "rand" trait

        # load trait vector (and other data we may or may not need later)
        load( 'simtrait.RData' )
        # loads: "causal_coeffs"  "causal_indexes" "trait"
        # calculations
        rmsds_rand[ rep ] <- trait_normality_rmsd( trait )
        mafs_rand[ rep ] <- causal_maf_min( X, causal_indexes )
        
        # process "inv" trait
        setwd( 'const_herit_loci' )
        
        # load trait vector (and other data we may or may not need later)
        load( 'simtrait.RData' )
        # loads: "causal_coeffs"  "causal_indexes" "trait"
        # calculations
        rmsds_inv[ rep ] <- trait_normality_rmsd( trait )
        mafs_inv[ rep ] <- causal_maf_min( X, causal_indexes )
        
        # go back down for next rep
        setwd( '../..' )
    }

    # return desired data
    return(
        list(
            rmsds_rand = rmsds_rand,
            rmsds_inv = rmsds_inv,
            mafs_rand = mafs_rand,
            mafs_inv = mafs_inv
        )
    )
}

############
### DATA ###
############

# move to where the data is
setwd( '../data/' )

# read dataset paths and names for output
datasets <- read_tsv( 'datasets.txt', col_types = 'cccii' )

# store data in a list, for boxplots
dataset_to_rmsds <- list()
dataset_to_mafs <- list()

for ( i in 1 : nrow( datasets ) ) {
    name_dir <- datasets$name_dir[ i ]
    setwd( name_dir )
    name_paper <- datasets$name_paper[ i ]
    message( name_paper)

    # decide if it's simulation or not, which decides where some genotype data ought to be
    is_sim <- grepl( 'sim', name_dir )

    # do all processing (combined rand/inv and rmsd/maf)
    obj <- process_reps( is_sim )
    
    # store "rand" traits
    name_out <- paste0( name_paper, ' RC' )
    dataset_to_rmsds[[ name_out ]] <- obj$rmsds_rand
    dataset_to_mafs[[ name_out ]] <- obj$mafs_rand

    # store "inv" traits
    name_out <- paste0( name_paper, ' FES' )
    dataset_to_rmsds[[ name_out ]] <- obj$rmsds_inv
    dataset_to_mafs[[ name_out ]] <- obj$mafs_inv
    
    # go back down
    setwd( '..' )
}

# visualize boxplots!
fig_start(
    'trait-normality',
    width = 6,
    height = 6,
    mar_b = 12
)
boxplot(
    dataset_to_rmsds,
    las = 2,
    ylab = 'RMSD( trait, norm-exp )'
)
fig_end()

# visualize boxplots!
fig_start(
    'causal-maf-mins',
    width = 6,
    height = 6,
    mar_b = 12
)
boxplot(
    dataset_to_mafs,
    las = 2,
    ylab = 'min causal MAF',
    log = 'y'
)
fig_end()
