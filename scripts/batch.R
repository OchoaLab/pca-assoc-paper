library(ochoalabtools)

# slurm job submission script for DCC

# shared items for all HGDP runs
#script <- 'real-05-gcta.R'
#plink <- FALSE
script <- 'real-06-pca-plink.R'
plink <- TRUE
bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3'
mem <- '4G' # probably a lot lower for PCA, but meh

# main submission steps
# global vars: script, bfile, mem, plink
submit_rep_pcs <- function(
                           rep,
                           pcs
                           ) {
    # pass params to command
    commands <- paste0(
        'time Rscript ', script,
        ' --bfile ', bfile,
        ' -r ', rep,
        ' --n_pcs ', pcs,
        ' --dcc'
    )
    
    # load plink module if needed
    if ( plink )
        commands <- c(
            'module load Plink/2.00a2LM',
            commands,
            'module unload Plink/2.00a2LM'
        )
    
    # give name clear params too
    name <- paste0( 'r', rep, 'p', pcs )
    
    # create submission file
    batch_writer(
        commands,
        name,
        mem = mem
    )

    # submit job!
    batch_submit( name )

    # remove script when done (not needed after submission)
    batch_cleanup( name )
}

# setup loop
reps_max <- 50
pcs_max <- 90

# I
# test one
submit_rep_pcs( rep = 1, pcs = 90 )

## # II
## # finish rest of run for these PCs
## pcs <- 90
## for ( rep in 2 : reps_max ) {
##     submit_rep_pcs( rep, pcs )
## }

## # III
## # loop (PCs on outside, as I was running plink locally)
## # only do things I was missing too
## pcs_list <- 57:89
## for ( pcs in pcs_list ) {
##     for ( rep in 1 : reps_max ) {
##         submit_rep_pcs( rep, pcs )
##     }
## }

## # IV
## # still missing a few for pc 56...
## pcs <- 56
## reps <- 35 : reps_max
## for ( rep in reps ) {
##     submit_rep_pcs( rep, pcs )
## }

########################

## # start loop (PCs on outside, as I was running plink locally)
## for ( pcs in 0 : pcs_max ) {
##     for ( rep in 1 : reps_max ) {
##         submit_rep_pcs( rep, pcs )
##     }
## }

## # start loop
## for ( rep in 1 : reps_max ) {
##     for ( pcs in 0 : pcs_max ) {
##         submit_rep_pcs( rep, pcs )
##     }
## }

