library(ochoalabtools)

# shared items for all HGDP runs
script <- 'real-05-gcta.R'
bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3'
mem <- '4G'
time <- NA

# setup loop
reps <- 50
pcs_max <- 90

# start loop
for ( rep in 1 : reps ) {
    for ( pcs in 0 : pcs_max ) {
        # pass params to command
        commands <- paste0(
            'time Rscript ', script, ' --bfile ', bfile, ' -r ', rep, ' --n_pcs ', pcs
        )
        # give name clear params too
        name <- paste0( 'hgdp_gcta_r', rep, '_pc', pcs )
        
        # create submission file
        batch_writer(
            commands,
            name,
            mem = mem,
            time = time
        )

        # submit job!
        batch_submit(
            name
        )
    }
}
