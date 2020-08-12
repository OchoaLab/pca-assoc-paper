library(ochoalabtools)

# shared items for all HGDP runs
script <- 'real-05-gcta.R'
bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3'
mem <- '4G'

# main submission steps
# global vars: script, bfile, mem
submit_rep_pcs <- function(
                           rep,
                           pcs
                           ) {
    # pass params to command
    commands <- paste0(
        'time Rscript ', script, ' --bfile ', bfile, ' -r ', rep, ' --n_pcs ', pcs
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
reps <- 50
pcs_max <- 90

# start loop
for ( rep in 1 : reps ) {
    for ( pcs in 0 : pcs_max ) {
        submit_rep_pcs( rep, pcs )
    }
}
