library(ochoalabtools)

# slurm job submission script for DCC

# shared items for all HGDP runs
plink <- FALSE
## plink <- TRUE
#bfile <- 'sim-n1000-k10-f0.1-s0.5-g1'; short <- 'l'
#bfile <- 'sim-n1000-k10-f0.1-s0.5-g20'; short <- 'f'
bfile <- 'HoPacAll_ld_prune_1000kb_0.3'; short <- 'h'
#bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3'; short <- 'd'
#bfile <- 'all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1'; short <- 'k'

# GCTA uses more memory, this works for the largest cases
mem <- if ( plink ) '4G' else '16G'

# select script automatically from boolean
script <- if ( plink ) 'real-06-pca-plink.R' else 'real-05-gcta.R'

# so names don't overlap between plink and gcta runs
short <- paste0( short, if (plink) 'p' else 'g' )

# main submission steps
# global vars: script, bfile, short, mem, plink
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
#        ' --sim'
    )
    
    # load plink module if needed
    if ( plink )
        # first add --plink flag (same line)
        commands <- paste0( commands, ' --plink' )
        # then load plink module (separate lines)
        commands <- c(
            'module load Plink/2.00a2LM',
            commands,
            'module unload Plink/2.00a2LM'
        )
    
    # give name clear params too
    name <- paste0( short, 'r', rep, 'p', pcs )
    
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
rep <- 1
pcs <- pcs_max
submit_rep_pcs( rep, pcs )

# II
# finish rest of rep
rep <- 1
for ( pcs in 0 : ( pcs_max - 1 ) ) {
    submit_rep_pcs( rep, pcs )
}

# III
# finish rest of reps
for ( rep in 2 : reps_max ) {
    for ( pcs in 0 : pcs_max ) {
        submit_rep_pcs( rep, pcs )
    }
}

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

