# gathers MAF distributions, for a plot of interest

library(BEDMatrix)
library(simtrait) # for allele_freqs

# constants
file_bed <- 'data'
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

# move to where the data is
setwd( '../data/' )

# get all the MAF vectors together in a list
mafs <- list()

# load each dataset
for ( dataset in datasets ) {
    # enter dir
    setwd( dataset )

    # slightly different processing for simulated data...
    is_sim <- grepl( 'sim-', dataset )
    
    # simulated datasets don't have a single genotype file, but one per rep
    # in that case enter rep-1
    if ( is_sim )
        setwd( 'rep-1' )
    
    # load genotypes BEDMatrix object
    X <- BEDMatrix( file_bed )

    # and get back down when done
    if ( is_sim )
        setwd( '..' )

    # calculate desired allele frequencies
    mafs[[ dataset ]] <- allele_freqs( X, fold = TRUE )

    # go back down for next dataset
    setwd( '..' )
}

# save all data on `data` base directory
save( mafs, file = 'mafs.RData' )
