library(ochoalabtools)

# slurm job submission script for DCC

# shared items for all HGDP runs
script <- 'real-05-gcta.R'
plink <- FALSE
#script <- 'real-06-pca-plink.R'
#plink <- TRUE
#bfile <- 'sim-n1000-k10-f0.1-s0.5-g1'; short <- 'l'
#bfile <- 'sim-n1000-k10-f0.1-s0.5-g20'; short <- 'f'
#bfile <- 'HoPacAll_ld_prune_1000kb_0.3'; short <- 'h'
#bfile <- 'hgdp_wgs_autosomes_ld_prune_1000kb_0.3'; short <- 'd'
bfile <- 'all_phase3_filt-minimal_ld_prune_1000kb_0.3_thinned-0.1'; short <- 'k'
mem <- '32G' # needed to increase for TGP GCTA
threads <- 4 # needed to make sure each process gets enough memory if all my jobs saturate the machines

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
        ' -t ', threads,
        ' --dcc' #,
#        ' --sim',
#        ' --plink'
    )
    
    # load plink module if needed
    if ( plink )
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
        mem = mem,
   	threads = threads
    )

    # submit job!
    batch_submit( name )

    # remove script when done (not needed after submission)
    batch_cleanup( name )
}

# setup loop
reps_max <- 50
pcs_max <- 90

# # I
# # test one
# submit_rep_pcs( rep = 1, pcs = 90 )
# submit_rep_pcs( rep = 1, pcs = 21 )
# submit_rep_pcs( rep = 1, pcs = 22 )

rep <- 1
# pcs_list <- 23:27
# pcs_list <- c( 33, 35, 37:39 )
pcs_list <- c( 40, 41, 43, 45:49 )
for ( pcs in pcs_list ) {
   submit_rep_pcs( rep = rep, pcs = pcs )
}

# # II
# # finish rest of run for these PCs
# pcs <- 90
# for ( rep in 2 : reps_max ) {
#     submit_rep_pcs( rep, pcs )
# }

# # III
# # loop (PCs on outside, as I was running plink locally)
# pcs_list <- 0:89
# for ( rep in 1 : reps_max ) {
#     for ( pcs in pcs_list ) {
#         submit_rep_pcs( rep, pcs )
#     }
# }

########################

## # start loop (PCs on outside, as I was running plink locally)
## for ( pcs in 0 : pcs_max ) {
##     for ( rep in 1 : reps_max ) {
##         submit_rep_pcs( rep, pcs )
##     }
## }

# # start loop
# for ( rep in 1 : reps_max ) {
#     for ( pcs in 0 : pcs_max ) {
#         submit_rep_pcs( rep, pcs )
#     }
# }

